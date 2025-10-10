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
}
bool is_stable(double dt, double dx, double dy){
    double mu_x = dt/(dx*dx);
    double mu_y = dt/(dy*dy);
    return mu_x + mu_y <= 0.5; // 稳定性条件
}
// 前向欧拉法求解二维Allen-Cahn热方程：u_t = \Delta u - f(u)/\epsilon^2.区域为[0,1]^2 * [0,T];周期边界条件
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
        vector<vector<vector<double>>> u;//依次为x，y，t
        vector<vector<vector<double>>> real_u;

    public:
        Allen_CahnEqSolver(double total_time, int grid_x, int grid_y, int time_steps, double epsilon):T(total_time),Nx(grid_x),Ny(grid_y),Nt(time_steps),ep(epsilon){
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
            double mu_t=dt/(ep*ep);
            for(int k=0;k<Nt;k++){
                for(int i=0;i<=Nx;i++){
                    for(int j=0;j<=Ny;j++){
                        int r = (i+1)%Nx;
                        int l = (i-1+Nx)%Nx;
                        int t = (j+1)%Ny;
                        int b = (j-1+Ny)%Ny;
                        u[i][j][k+1] = (1.0-2*(mu_x+mu_y))*u[i][j][k]+mu_x*(u[r][j][k]+u[l][j][k])+mu_y*(u[i][t][k]+u[i][b][k])
                        -mu_t*source(u[i][j][k]);
                    }
                }
            }
        }
        
        
        const vector<vector<vector<double>>>& getsolution() const{
            return u;
        }

};
int main(){
    double T = 0.1;
    int Nx = 10;
    int Ny = 10;
    int Nt = 1000;
    double ep = 0.1;
    if(!is_stable(T/Nt, 1.0/Nx, 1.0/Ny)){
        cout << "The method is not stable with the given parameters." << endl;
        system("pause");
        return 0;
    }
    Allen_CahnEqSolver solver(T, Nx, Ny, Nt, ep);
    solver.solve();
    auto u = solver.getsolution();
    // 输出结果
    for(int k=0;k<=Nt;k++){
        for(int i=0;i<=Nx;i++){
            for(int j=0;j<=Ny;j++){
                cout << "u(" << i*1.0/Nx << "," << j*1.0/Ny << "," << T*k*1.0/Nt << ") = " << u[i][j][k] << endl;
            }
        }
    }
    system("pause");
    return 0;
}
