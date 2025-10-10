#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <complex>
#include <new>  // 对于 ::operator new
#include <memory> // 对于标准分配器接口
// 周期边界条件下求解 -\Delta u + u = f
using namespace std;
using Func = function<double(double,double)>;
const double PI = 3.1415926535897932384626;
Func f = [](double x,double y){
        return (8*PI*PI+1)*sin(2*PI*x)*sin(2*PI*y);
};
Func g = [](double x,double y){
        return 0.0;
};
Func exact_solution = [](double x,double y){
        return sin(2*PI*x)*sin(2*PI*y);
};
class PoissonSolver{
    private:
        int N;//剖分数，格点数=N+1
        double h=1.0/N;
        Func f;
        Func g;
        vector<vector<double>> u;
        
        double dotProduct (const vector<double>& a, const vector<double>& b){
            double result = 0.0;
            for(int i = 0;i < a.size();i++){
                result += a[i]*b[i];
            }
            return result;
        }
        // 用FFT求解线性方程组
        vector<complex<double>> FFT(const vector<complex<double>>& x){
            int n = x.size();//设定n为2的幂次方
            // FFT算法实现
            if (n == 1) return x;
            vector<complex<double>> x_even(n/2), x_odd(n/2);
            for (int i = 0; i < n/2; i++) {
                x_even[i] = x[2*i];
                x_odd[i] = x[2*i + 1];
            }
            vector<complex<double>> fx_even = FFT(x_even);
            vector<complex<double>> fx_odd = FFT(x_odd);
            vector<complex<double>> fx(n);
            fill(fx.begin(), fx.end(), complex<double>(0.0, 0.0));
            vector<complex<double>> w(n/2);
            for(int k=0;k<n/2;k++)
                w[k] = exp(complex<double>(0, -2*PI*k/n));
            for (int k = 0; k < n/2; k++) {
                fx[k] = fx_even[k] + w[k] * fx_odd[k];
                fx[k + n/2] = fx_even[k] - w[k] * fx_odd[k];
            }
            return fx;  
        }
        vector<complex<double>> IFFT(const vector<complex<double>>& fx){
            int n = fx.size();
            // IFFT算法实现
            if (n == 1) return fx;
            vector<complex<double>> fx_even(n/2), fx_odd(n/2);
            for (int i = 0; i < n/2; i++) {
                fx_even[i] = fx[2*i];
                fx_odd[i] = fx[2*i + 1];
            }
            vector<complex<double>> x_even = IFFT(fx_even);
            vector<complex<double>> x_odd = IFFT(fx_odd);
            vector<complex<double>> x(n);
            fill(x.begin(), x.end(), complex<double>(0.0, 0.0));
            vector<complex<double>> w(n/2);
            for(int k=0;k<n/2;k++)
                w[k] = exp(complex<double>(0, 2*PI*k/n));
            for (int k = 0; k < n/2; k++) {
                x[k] = (x_even[k] + w[k] * x_odd[k])/double(n);
                x[k + n/2] = (x_even[k] - w[k] * x_odd[k])/double(n);
            }
            return x;  
        }
        // 2D FFT实现
    vector<vector<complex<double>>> FFT2D(const vector<vector<complex<double>>>& input) {
        int rows = input.size();
        int cols = input[0].size();
        vector<vector<complex<double>>> output(rows, vector<complex<double>>(cols));
        
        // 对每一行做FFT
        vector<vector<complex<double>>> row_fft(rows, vector<complex<double>>(cols));
        for (int i = 0; i < rows; i++) {
            row_fft[i] = FFT(input[i]);
        }
        
        // 转置矩阵
        vector<vector<complex<double>>> transposed(cols, vector<complex<double>>(rows));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = row_fft[i][j];
            }
        }
        
        // 对每一列做FFT
        for (int j = 0; j < cols; j++) {
            auto col_fft = FFT(transposed[j]);
            for (int i = 0; i < rows; i++) {
                output[i][j] = col_fft[i];
            }
        }
        
        return output;
    }
    
    // 2D IFFT实现
    vector<vector<complex<double>>> IFFT2D(const vector<vector<complex<double>>>& input) {
        int rows = input.size();
        int cols = input[0].size();
        vector<vector<complex<double>>> output(rows, vector<complex<double>>(cols));
        
        // 对每一行做IFFT
        vector<vector<complex<double>>> row_ifft(rows, vector<complex<double>>(cols));
        for (int i = 0; i < rows; i++) {
            row_ifft[i] = IFFT(input[i]);
        }
        
        // 转置矩阵
        vector<vector<complex<double>>> transposed(cols, vector<complex<double>>(rows));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = row_ifft[i][j];
            }
        }
        
        // 对每一列做IFFT
        for (int j = 0; j < cols; j++) {
            auto col_ifft = IFFT(transposed[j]);
            for (int i = 0; i < rows; i++) {
                output[i][j] = col_ifft[i];
            }
        }
        
        return output;
    }
    public:
        
        PoissonSolver(int grid_size, Func source, Func boundary):N(grid_size),f(source),g(boundary){
            h = 1.0/N;
            u.resize(N,vector<double>(N,0.0));
        } 
        void solve(){
        // 创建右端项矩阵
        vector<vector<complex<double>>> B(N, vector<complex<double>>(N, complex<double>(0.0, 0.0)));
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double x =  h*i;  // 网格点位置
                double y =  h*j;
                B[i][j] = complex<double>(h*h*f(x, y), 0.0);
            }
        }
        
        auto B_hat = FFT2D(B);
        vector<vector<complex<double>>> Lambda(N, vector<complex<double>>(N, complex<double>(0.0, 0.0)));
        vector<vector<complex<double>>> U_hat(N, vector<complex<double>>(N, complex<double>(0.0, 0.0)));
        for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
                int k1 = (k < N/2) ? k : k - N; // 处理周期性边界条件
                int l1 = (l < N/2) ? l : l - N; // 处理周期性边界条件
                complex<double> lambda = complex<double> (4.0 + h*h - 2*cos(2*PI*k1/N) - 2*cos(2*PI*l1/N),0.0); // 特征值
                Lambda[k][l] = lambda;
                U_hat[k][l] = B_hat[k][l] / Lambda[k][l];

            }
        }
        
        
        auto U = IFFT2D(U_hat);
        vector<vector<double>> u(N, vector<double>(N));
        // 存储结果
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                u[i][j] = U[i][j].real();
            }
        }
        }
        const vector<vector<double>>& getSolution() const {
            return u;
        }
};
int main(){
    int N = 8;
    double h = 1.0/N;
    PoissonSolver solver(N,f,g);
    solver.solve();
    auto u = solver.getSolution();
    vector<double> x_grid(N),y_grid(N);
    for(int i=0;i<N;i++){
        x_grid[i] = i*h;
        y_grid[i] = i*h;
    }
    vector<vector<double>> exact_u,error;
    exact_u.resize(N,vector<double>(N,0.0));
    error.resize(N,vector<double>(N,0.0));
    double max_error = 0.0;
    for (int i=0;i<N;i++){ 
        for (int j=0;j<N;j++){
            exact_u[i][j] = exact_solution(x_grid[i], y_grid[j]);
            error[i][j] = abs(u[i][j] - exact_u[i][j]);
            cout << error[i][j] << " ";
            max_error = max(max_error, error[i][j]);
        }
    }
    cout << "Max error: " << max_error << endl;
    system("pause");
    return 0;
}
