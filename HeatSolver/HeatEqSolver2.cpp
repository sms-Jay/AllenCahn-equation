#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
using namespace std;
using Func = function<double(double,double)>;
const double PI = 3.1415926535897932384626;
Func initial = [](double x,double y){
    return (1+sin(2*PI*x)*sin(2*PI*y))/3.0;
};
double solu(double x,double y,double t){
    return (1+sin(2*PI*x)*sin(2*PI*y)*exp(-8*PI*PI*t))/3.0;
}
// 后向欧拉法求解二维热方程：u_t = \Delta u.区域为[0,1]^2 * [0,T];周期边界条件
// 我们采取共轭梯度法；
class HeatEqSolver{
    private:
        double T;
        int Nx;
        int Ny;
        int Nt;
        double dx;
        double dy;
        double dt;
        vector<vector<vector<double>>> u;//依次为t,x，y
        vector<vector<vector<double>>> real_u;

    public:
        HeatEqSolver(double total_time, int grid_x, int grid_y, int time_steps):T(total_time),Nx(grid_x),Ny(grid_y),Nt(time_steps){
            dx = 1.0/Nx;
            dy = 1.0/Ny;
            dt = T/Nt;
            u.resize(Nt+1,vector<vector<double>>(Nx,vector<double>(Ny,0.0)));
            real_u.resize(Nt+1, vector<vector<double>>(Nx, vector<double>(Ny, 0.0)));
            // 初始化初始条件
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    double x = i*dx;
                    double y = j*dy;
                    u[0][i][j] = initial(x,y);
                }
            }
        }

        double dot_product(const vector<double>& a, const vector<double>& b) {
            double result = 0.0;
            int n = a.size();
            for (size_t i = 0; i < n; i++) {
                result += a[i] * b[i];
            }
            return result;
        }

        vector<double> matrix_vector_product(const vector<double>& x){
            vector<double> Ax(x.size(), 0.0);
            double mu_x=dt/(dx*dx);
            double mu_y=dt/(dy*dy);
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    int idx = i*Ny+j;
                    Ax[idx] = (1.0 + 2*(mu_x + mu_y))*x[idx];
                    int r = (i+1)%Nx;
                    int l = (i-1+Nx)%Nx;
                    int t = (j+1)%Ny;
                    int b = (j-1+Ny)%Ny;
                    Ax[idx] -= mu_x*(x[r*Ny+j] + x[l*Ny+j]);
                    Ax[idx] -= mu_y*(x[i*Ny+t] + x[i*Ny+b]);
                }
            }
            return Ax;
        }
        
        void CG(const vector<double>& b, vector<double>& x, double tol, int max_iter){
            int n = b.size();
            vector<double> r = b;
            vector<double> p = r;
            vector<double> Ap(n);
            double b_norm = sqrt(dot_product(b,b));
            if (b_norm < 1e-10) b_norm = 1.0;
            double rold = dot_product(r, r);
            
            for (int k = 0; k < max_iter;k++){
                auto Ap = matrix_vector_product(p);
                double pAp = dot_product(p,Ap);
                double alpha = rold/pAp;
                for (int i=0;i<n;i++){
                    x[i] +=alpha*p[i];
                    r[i] -=alpha*Ap[i];
                }
                double rnew = dot_product(r,r);
                double r_norm = sqrt(rnew);
                if(r_norm/b_norm < tol) break;
                double beta = rnew/rold;
                for (int i=0;i<n;i++){
                    p[i] = r[i] + beta*p[i];
                }
                rold = rnew;
            }
        }
        
        vector<double> RHS(vector<vector<double>>& x){
            int n = Nx*Ny;
            vector<double> b(n, 0.0);
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    int idx = i*Ny+j;
                    b[idx] = x[i][j];
                }
            }
            return b;
        }

        void solve(){
            int n = Nx*Ny;
            double mu_x=dt/(dx*dx);
            double mu_y=dt/(dy*dy);
            for(int k = 1;k<=Nt;k++){
                vector<double> b(n,0.0),x(n,0.0);
                b = RHS(u[k-1]);
                CG(b,x,1e-8,1000000);
                for(int i=0;i<Nx;i++){
                    for(int j=0;j<Ny;j++){
                        int idx = i*Ny+j;
                        u[k][i][j]=x[idx];
                    }
                }
            }
        }
        
        void realsolu(){
            for(int k=0;k<=Nt;k++){
                for(int i=0;i<Nx;i++){
                    for(int j=0;j<Ny;j++){
                        double x = i*dx;
                        double y = j*dy;
                        double t = k*dt;
                        real_u[k][i][j] = solu(x,y,t);
                    }
                }
            }
        }
        const vector<vector<vector<double>>>& getsolution() const{
            return u;
        }
        const vector<vector<vector<double>>>& getrealsolution() const{
            return real_u;
        }
};
int main(){
    double T = 0.1;
    int Nx = 100;
    int Ny = 100;
    int Nt = 10000;
    HeatEqSolver solver(T, Nx, Ny, Nt);
    solver.solve();
    solver.realsolu();
    auto u = solver.getsolution();
    auto real_u = solver.getrealsolution();
    double max_error = 0.0;
    for(int k=0;k<=Nt;k++){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                double error = abs(u[k][i][j] - real_u[k][i][j]);
                max_error = max(max_error, error);
            }
        }
    }
    cout << "Max error: " << max_error << endl;
    system("pause");
    return 0;
}
