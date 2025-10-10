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
double source(double u){
    return log(1+u)-log(1-u)-3*u;
};
// 向后欧拉法求解二维Allen-Cahn热方程：u_t = \Delta u - f(u)/\epsilon^2.区域为[0,1]^2 * [0,T];周期边界条件
class Allen_CahnEqSolver{
    private:
        double T;
        int Nx;
        int Ny;
        int Nt;
        double dx;
        double dy;
        double dt;
        double ep; // 控制非线性项的强度
        vector<vector<vector<double>>> u;//依次为t,x,y

    public:
        Allen_CahnEqSolver(double total_time, int grid_x, int grid_y, int time_steps, double epsilon):T(total_time),Nx(grid_x),Ny(grid_y),Nt(time_steps),ep(epsilon){
            dx = 1.0/Nx;
            dy = 1.0/Ny;
            dt = T/Nt;
            u.resize(Nt+1,vector<vector<double>>(Nx,vector<double>(Ny,0.0)));
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
        vector<double> elementwise_mult(const vector<double>& a,const vector<double>& b){
            int n = a.size();
            vector<double> c(n,0.0);
            for(int i=0;i<n;i++){
                c[i]=a[i]*b[i];
            }
            return c;
        }
        vector<double> vec_add(const vector<double>&a, const vector<double>& b){
            int n=a.size();
            vector<double> c(n,0.0);
            for(int i=0;i<n;i++){
                c[i]=a[i]+b[i];
            }
            return c;
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
        bool is_in(const vector<double>& x){
            int n = x.size();
            for(int i=0;i<n;i++){
                if(x[i]>=1.0 || x[i]<=-1.0) return false;
            }
            return true;
        }
        vector<double> vec_mul(const double a, const vector<double>& x){
            int n = x.size();
            vector<double> y(n,0.0);
            for(int i=0;i<n;i++){
                y[i]=1.0*a*x[i];
            }
            return y;
        }
        vector<double> gradient_f(const vector<double>& x) {
            int n = x.size();
            double mu_t=dt/(ep*ep);
            vector<double> gf(n, 0.0);
            for (int i = 0; i < n; i++) {
                gf[i] = mu_t*2.0/(1-x[i]*x[i]) - 3.0;
            }
            return gf;
        }
        void CG(const vector<double>& b, vector<double>& x, vector<double>& dx, double tol, int max_iter){
            int n = b.size();
            vector<double> r = b;
            vector<double> p = r;
            vector<double> Ap(n,0.0);
            vector<double> gf = gradient_f(x);
            double b_norm = sqrt(dot_product(b,b));
            if (b_norm < 1e-10) b_norm = 1.0;
            double rold = dot_product(r, r);

            for (int k = 0; k < max_iter;k++){
                auto Ap = matrix_vector_product(p);
                auto c = elementwise_mult(gf,p);
                Ap = vec_add(Ap,c);
                double pAp = dot_product(p,Ap);
                double alpha = rold/pAp;
                for (int i=0;i<n;i++){
                    dx[i] +=alpha*p[i];
                    r[i] -=alpha*Ap[i];
                }
                double rnew = dot_product(r,r);
                double r_norm = sqrt(rnew);
                if(r_norm/b_norm < tol){
                    // cout<<"Converged in "<<k+1<<" iterations."<<endl;
                    break;
                }
                double beta = rnew/rold;
                for (int i=0;i<n;i++){
                    p[i] = r[i] + beta*p[i];
                }
                rold = rnew;
                if(k == max_iter-1){
                    cout<<"Warning: CG did not converge within the maximum iterations."<<endl;
                }
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
        
        vector<double> nonlinear(const vector<double>& x){
            //计算source作用下的结果
            int n = Nx*Ny;
            double mu_t=dt/(ep*ep);
            vector<double> f(n,0.0);
            for(int i=0;i<n;i++){
                f[i] = mu_t*source(x[i]);
            }
            return f;
        }
        void newton(vector<double>& x,const vector<double>& b, int max_iter, double tol) {
        //牛顿法求解方程Ax+f(x)=b，初值选取为上一时刻的终值;F(x) = Ax + f(x) - b
            int n=Nx*Ny;
            for (int iter = 0; iter < max_iter; iter++) {
                vector<double> Ax = matrix_vector_product(x);
                vector<double> f = nonlinear(x);
                // 计算残差
                vector<double> r(n, 0.0);
                for (int i = 0; i < n; i++) {
                    r[i] = -Ax[i] - f[i] + b[i];
                }
                
                // 检查收敛性
                double norm_r = sqrt(dot_product(r, r));
                if (norm_r < tol) {
                    break;
                }
                vector<double> dx(n,0.0);
                //共轭梯度法求解：(gradient_f(x) + A)*dx = r
                CG(r,x,dx,1e-4,10000);
                /*
                for(double i = 1.0;i>=-1.0;i-=0.01){
                    auto y = vec_add(x,vec_mul(i,x));
                    if(is_in(y)){
                        vector<double> Ay = matrix_vector_product(y);
                        vector<double> fy = nonlinear(y);
                        // 计算残差
                        vector<double> ry(n, 0.0);
                        for (int i = 0; i < n; i++) {
                            ry[i] = b[i] - Ay[i] - fy[i];
                        }
                        double norm_ry = sqrt(dot_product(ry, ry));
                        if(norm_ry < norm_r/2){
                            cout<<i<<endl;
                            x=y;
                            break;
                        }
                    }
                }
                */
                // 确保解在(-1,1)范围内
                
                x=vec_add(x,vec_mul(0.1,dx));
                /*
                for (int i = 0; i < x.size(); i++) {
                    if (x[i] >= 1.0) x[i] =0.9999;
                    if(x[i]<=-1.0) x[i]=-0.9999;
                }
                */
                if(iter == max_iter-1){
                    cout<<"Warning: Newton's method did not converge within the maximum iterations."<<endl;
                }
            }   
            
        }
        void solve(){
            double mu_x=dt/(dx*dx);
            double mu_y=dt/(dy*dy);
            for(int k=1;k<=Nt;k++){
                vector<double> b = RHS(u[k-1]);
                vector<double> x = b;//本迭代初值为上一时刻终值
                newton(x,b,10000,1e-4);
                for(int i=0;i<Nx;i++){
                    for(int j=0;j<Ny;j++){
                        int idx=i*Ny+j;
                        u[k][i][j]=x[idx];
                    }
                }
            }
        }
        
        const vector<vector<vector<double>>>& getsolution() const{
            return u;
        }

};
int main(){
    double T = 1e10;
    int Nx = 20;
    int Ny = 20;
    int Nt = 10;
    double ep = 0.05;
    Allen_CahnEqSolver solver(T, Nx, Ny, Nt, ep);
    solver.solve();
    auto u = solver.getsolution();
    // 输出结果
    for(int k=Nt;k<=Nt;k++){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                cout << "u(" << T*k*1.0/Nt << "," << i*1.0/Nx << "," << j*1.0/Ny << ") = " << u[k][i][j] << endl;
            }
        }
    }
    system("pause");
    return 0;
}
