#ifndef FFT_2D_SOLVER_H
#define FFT_2D_SOLVER_H

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <omp.h>

class FFT2DSolver {
private:
    int Nx, Ny, N;
    double dx, dy, dx2, dy2;
    std::vector<std::complex<double>> cached_laplace_eigenvalues;
    bool eigenvalues_cached;

    const double PI = 3.1415926535897932384626;

public:
    // 构造函数
    FFT2DSolver(int grid_x, int grid_y, double delta_x, double delta_y) 
        : Nx(grid_x), Ny(grid_y), dx(delta_x), dy(delta_y) {
        N = Nx * Ny;
        dx2 = dx * dx;
        dy2 = dy * dy;
        eigenvalues_cached = false;
    }

    // 预计算拉普拉斯特征值（周期边界）
    void precompute_laplace_eigenvalues() {
        if (eigenvalues_cached) return;
        
        cached_laplace_eigenvalues.resize(N);
        
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < Nx; i++) {
            for(int j = 0; j < Ny; j++) {
                int idx = i * Ny + j;
                
                int kx = (i <= Nx/2) ? i : i - Nx;
                int ky = (j <= Ny/2) ? j : j - Ny;
                
                cached_laplace_eigenvalues[idx] = std::complex<double>(
                    -2.0 * (
                        (1.0 - std::cos(2.0 * PI * kx / Nx)) / dx2 + 
                        (1.0 - std::cos(2.0 * PI * ky / Ny)) / dy2
                    ), 0.0
                );
            }
        }
        
        eigenvalues_cached = true;
    }

    // 获取拉普拉斯特征值
    const std::vector<std::complex<double>>& get_laplace_eigenvalues() {
        precompute_laplace_eigenvalues();
        return cached_laplace_eigenvalues;
    }

    // 1D FFT (Cooley-Tukey算法)
    void fft_1d(std::vector<std::complex<double>>& x, bool inverse) const {
        int n = x.size();
        if (n <= 1) return;
        
        // 位反转重排
        for (int i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; j >= bit; bit >>= 1) {
                j -= bit;
            }
            j += bit;
            if (i < j) {
                std::swap(x[i], x[j]);
            }
        }
        
        // 迭代计算FFT
        for (int len = 2; len <= n; len <<= 1) {
            double angle = 2 * PI / len * (inverse ? 1 : -1);
            std::complex<double> wlen(std::cos(angle), std::sin(angle));
            
            for (int i = 0; i < n; i += len) {
                std::complex<double> w(1);
                for (int j = 0; j < len/2; j++) {
                    std::complex<double> u = x[i + j];
                    std::complex<double> v = w * x[i + j + len/2];
                    x[i + j] = u + v;
                    x[i + j + len/2] = u - v;
                    w *= wlen;
                }
            }
        }
        
        // 逆变换需要归一化
        if (inverse) {
            for (int i = 0; i < n; i++) {
                x[i] /= n;
            }
        }
    }

    // 2D FFT
    std::vector<std::complex<double>> fft_2d(const std::vector<std::complex<double>>& input, bool inverse) const {
        std::vector<std::complex<double>> output(N);
        std::vector<std::complex<double>> temp(N);
        
        // 对每一行进行FFT
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < Nx; i++) {
            std::vector<std::complex<double>> row(Ny);
            for (int j = 0; j < Ny; j++) {
                row[j] = input[i * Ny + j];
            }
            fft_1d(row, inverse);
            for (int j = 0; j < Ny; j++) {
                temp[i * Ny + j] = row[j];
            }
        }
        
        // 对每一列进行FFT
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < Ny; j++) {
            std::vector<std::complex<double>> col(Nx);
            for (int i = 0; i < Nx; i++) {
                col[i] = temp[i * Ny + j];
            }
            fft_1d(col, inverse);
            for (int i = 0; i < Nx; i++) {
                output[i * Ny + j] = col[i];
            }
        }
        
        return output;
    }

    // 实数输入的2D FFT（自动转换为复数）
    std::vector<std::complex<double>> fft_2d_real(const std::vector<double>& input, bool inverse) const {
        std::vector<std::complex<double>> complex_input(N);
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < N; i++) {
            complex_input[i] = std::complex<double>(input[i], 0.0);
        }
        return fft_2d(complex_input, inverse);
    }

    // 求解线性系统: (alpha * I - beta * Laplace) * x = b
    std::vector<double> solve_linear_system(const std::vector<double>& b, 
                                          double alpha, double beta = 1.0) {
        precompute_laplace_eigenvalues();
        
        std::vector<double> x(N, 0.0);
        
        // 转换到频率域
        auto b_freq = fft_2d_real(b, false);
        
        // 在频率域求解
        std::vector<std::complex<double>> x_freq(N);
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < N; i++) {
            double system_eigenvalue = alpha - beta * cached_laplace_eigenvalues[i].real();
            
            if(std::abs(system_eigenvalue) > 1e-12) {
                x_freq[i] = b_freq[i] / system_eigenvalue;
            } else {
                x_freq[i] = 0.0;
            }
        }
        
        // 转换回空间域
        auto x_complex = fft_2d(x_freq, true);
        
        // 提取实部
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < N; i++) {
            x[i] = x_complex[i].real();
        }
        
        return x;
    }

    // 获取网格信息
    int get_total_size() const { return N; }
    int get_nx() const { return Nx; }
    int get_ny() const { return Ny; }
    
    // 清理缓存
    void clear_cache() {
        cached_laplace_eigenvalues.clear();
        eigenvalues_cached = false;
    }
};

#endif // FFT_2D_SOLVER_H