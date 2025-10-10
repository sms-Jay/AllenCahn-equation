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
bool is_stable(double dt, double dx, double dy){
    double mu_x = dt/(dx*dx);
    double mu_y = dt/(dy*dy);
    return mu_x + mu_y <= 0.5; // 稳定性条件
}
// 后向欧拉法求解二维热方程：u_t = \Delta u.区域为[0,1]^2 * [0,T];周期边界条件
class HeatEqSolver{
    private:
        double T;
        int Nx;
        int Ny;
        int Nt;
        double dx;
        double dy;
        double dt;
        vector<vector<vector<double>>> u;//依次为x，y，t
        vector<vector<vector<double>>> real_u;

    public:
        HeatEqSolver(double total_time, int grid_x, int grid_y, int time_steps):T(total_time),Nx(grid_x),Ny(grid_y),Nt(time_steps){
            dx = 1.0/Nx;
            dy = 1.0/Ny;
            dt = T/Nt;
            u.resize(Nx+1,vector<vector<double>>(Ny+1,vector<double>(Nt+1,0.0)));
            real_u.resize(Nx+1, vector<vector<double>>(Ny+1, vector<double>(Nt+1, 0.0)));
            // 初始化初始条件
            for(int i=0;i<=Nx;i++){
                for(int j=0;j<=Ny;j++){
                    double x = i*dx;
                    double y = j*dy;
                    u[i][j][0] = initial(x,y);
                }
            }
        }
        void solve(){
            double mu_x=dt/(dx*dx);
            double mu_y=dt/(dy*dy);
            for(int k=0;k<Nt;k++){
                for(int i=0;i<=Nx;i++){
                    for(int j=0;j<=Ny;j++){
                        int r = (i+1)%Nx;
                        int l = (i-1+Nx)%Nx;
                        int t = (j+1)%Ny;
                        int b = (j-1+Ny)%Ny;
                        u[i][j][k+1] = (1.0-2*(mu_x+mu_y))*u[i][j][k]+mu_x*(u[r][j][k]+u[l][j][k])+mu_y*(u[i][t][k]+u[i][b][k]);
                    }
                }
            }
        }
        
        void realsolu(){
            for(int k=0;k<=Nt;k++){
                for(int i=0;i<=Nx;i++){
                    for(int j=0;j<=Ny;j++){
                        double x = i*dx;
                        double y = j*dy;
                        double t = k*dt;
                        real_u[i][j][k] = solu(x,y,t);
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
    if(!is_stable(T/Nt, 1.0/Nx, 1.0/Ny)){
        cout << "The method is not stable with the given parameters." << endl;
        system("pause");
        return 0;
    }
    HeatEqSolver solver(T, Nx, Ny, Nt);
    solver.solve();
    solver.realsolu();
    auto u = solver.getsolution();
    auto real_u = solver.getrealsolution();
    double max_error = 0.0;
    for(int i=0;i<=Nx;i++){
        for(int j=0;j<=Ny;j++){
            for(int k=0;k<=Nt;k++){
                double error = abs(u[i][j][k] - real_u[i][j][k]);
                max_error = max(max_error, error);
            }
        }
    }
    cout << "Max error: " << max_error << endl;
    system("pause");
    return 0;
}
