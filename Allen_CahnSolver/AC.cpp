#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
using namespace std;
const double PI = 3.1415926535897932384626;
double ini_1(double x, double y){
    return (1+sin(2*PI*x)*sin(2*PI*y))/3.0;
}
double ini_2(double x, double y){
    if(abs(x)<=0.35 && abs(y)<=0.35) return (1e-5);
    else return 1-(1e-5);
}
double ini_3(double x, double y){
    return (double) rand() / (RAND_MAX + 1.0);
}
double F(double t,double theta){
    if (t>0 && t<1e-16) return (1-t)*log(1-t)+theta*t*(1-t);
    else if (t<1 && t>1-1e-16) return t*log(t)+theta*t*(1-t);
    else if (t>1e-16 && t<1-1e-16) return t*log(t)+(1-t)*log(1-t)+theta*t*(1-t);
    else{
        cout<<"Wrong range!"<<endl;
        return 0.0;
    }
}
void saveDataToFile(const vector<vector<vector<double>>>& u, const string& filename) {
    ofstream outFile(filename);
    
    if (!outFile.is_open()) {
        cerr << "Cannot Open File: " << filename <<endl;
        return;
    }
    
    outFile << fixed <<setprecision(6);
    
    int timeSteps = u.size();
    int xSize = u[0].size();
    int ySize = u[0][0].size();
    
    outFile << timeSteps << " " << xSize << " " << ySize << endl;
    
    for (int x = 0; x < xSize; ++x) {
        outFile << static_cast<double>(x);
        if (x < xSize - 1) outFile << " ";
    }
    outFile << endl;
    
    for (int y = 0; y < ySize; ++y) {
        outFile << static_cast<double>(y);
        if (y < ySize - 1) outFile << " ";
    }
    outFile << endl;
    
    for (int t = 0; t < timeSteps; ++t) {
        outFile << "t=" << t << endl;
        for (int x = 0; x < xSize; ++x) {
            for (int y = 0; y < ySize; ++y) {
                outFile << u[t][x][y];
                if (y < ySize - 1) outFile << " ";
            }
            outFile << endl;
        }
    }
    outFile.close();
    cout << "Data Saved to: " << filename << endl;
}
class AC{
    private:
        int Nt;
        double dt;
        int Nx;
        double dx;
        double dx2;
        int Ny;
        double dy;
        double dy2;
        double dxdy;
        int N;
        double ep;
        double ep2;
        double theta=4.0;
        vector<vector<double>> u;
        vector<vector<vector<double>>> U;
        vector<double> energy;

    public:
        AC(double time, int time_steps, int grid_x, int grid_y,  double epsilon):dt(time),Nt(time_steps),Nx(grid_x),Ny(grid_y),ep(epsilon){
            dx = 2.0/Nx;
            dx2 = dx*dx;
            dy = 2.0/Ny;
            dy2 = dy*dy;
            dxdy=dx*dy;
            ep2 = ep*ep;
            N = Nx*Ny;
            u.resize(Nt+1,vector<double>(N,0.0));
            U.resize(Nt+1,vector<vector<double>>(Nx,vector<double>(Ny,0.0)));
            energy.resize(Nt+1,0.0);
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    double x = i*dx-1.0;
                    double y = j*dy-1.0;
                    int idx = i*Ny+j;
                    u[0][idx] = ini_3(x,y);
                }
            }
        }
        
        double dotProduct (const vector<double>& a, const vector<double>& b){
            double result = 0.0;
            for(int i = 0;i < a.size();i++){
                result += a[i]*b[i];
            }
            return result;
        }

        double Energy(const vector<double>& U){
            double E = 0.0;
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    int ip = (i+1)%Nx;
                    int im = (i-1+Nx)%Nx;
                    int jp = (j+1)%Ny;
                    int jm = (j-1+Ny)%Ny;
                    int idx = i*Ny+j;
                    int idx_ip = ip*Ny+j;
                    int idx_im = im*Ny+j;
                    int idx_jp = i*Ny+jp;
                    int idx_jm = i*Ny+jm;
                    double gx = (U[idx_ip]-U[idx_im])/(2*dx);// 
                    double gy = (U[idx_jp]-U[idx_jm])/(2*dy);
                    E += dxdy*(0.5*ep2*(gx*gx+gy*gy) + F(U[idx],theta));
                }
            }
            return E;
        }

        vector<double> matrixVecProduct(const vector<double>& U, double rho){
            // calculate£¨1/dt -ep^2 \Delta + rho)*U
            auto AU = U;
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    int idx = i*Ny+j;
                    AU[idx] = (1.0/dt + rho)*U[idx];
                    int ip = (i+1)%Nx;
                    int im = (i-1+Nx)%Nx;
                    int jp = (j+1)%Ny;
                    int jm = (j-1+Ny)%Ny;
                    int idx_ip = ip*Ny+j;
                    int idx_im = im*Ny+j;
                    int idx_jp = i*Ny+jp;
                    int idx_jm = i*Ny+jm;
                    AU[idx] -= ep2*( (U[idx_ip]-2*U[idx]+U[idx_im])/(dx2) + (U[idx_jp]-2*U[idx]+U[idx_jm])/(dy2) );
                }
            }
            return AU;
        }
        
        vector<double> rhs(const vector<double>& Un, const vector <double>& Y, double rho,const vector<double>& U2){
            auto b = Un;
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    int idx = i*Ny+j;
                    b[idx] = Un[idx]/dt -Y[idx] +rho*U2[idx];
                }
            }
            return b;
        }
        // Symmetric and Positive definenite, consider FFT or PCG
        vector<double> CG(const vector<double>& b, vector<double>& x, double tol, double rho, int max_iter){
            vector<double> r = b;
            vector<double> p = r;
            vector<double> Ap(N);
            double b_norm = sqrt(dotProduct(b,b));
            if (b_norm < 1e-10) b_norm = 1.0;
            double rold = dotProduct(r, r);
            
            for (int k = 1; k <= max_iter;k++){
                auto Ap = matrixVecProduct(p,rho);
                double pAp = dotProduct(p,Ap);
                double alpha = rold/pAp;
                for (int i=0;i<N;i++){
                    x[i] +=alpha*p[i];
                    r[i] -=alpha*Ap[i];
                }
                double rnew = dotProduct(r,r);
                double r_norm = sqrt(rnew);
                if(r_norm/b_norm < tol) {
                    // cout<<"CG converge at "<<k<<endl;
                    break;
                }
                double beta = rnew/rold;
                for (int i=0;i<N;i++){
                    p[i] = r[i] + beta*p[i];
                }
                rold = rnew;
                if(k==max_iter) cout<<"CG not converge."<<endl;
            }
            return x;
        }
        
        double newton(const double ini, const double u1, const double un, const double y, const double rho){
            double u = ini;
            for(int iter=1;iter<=100;iter++){
                double fu = log(u) - log(1-u) + theta*(1-2*un) - y - rho*(u1 - u);
                double dfu = 1.0/(u*(1-u)) + rho;
                double u_new = u - fu/dfu;
                if(fabs(f(u_new,u1,un,y,rho)) < 1e-8){
                    // cout<<"Newton converge at "<<iter<<endl;
                    return u_new;
                }
                if(iter == 100) {
                    cout<<"Newton not converge."<<endl;
                }
                u = u_new;
            }
            return u;
        }
         
        double f(const double u2, const double u1, const double un, const double y, const double rho){
            return log(u2/(1-u2)) + theta*(1-2*un) - y - rho*(u1 - u2);
        }

        double solveu2(const double u2, const double u1, const double un, const double y, const double rho){
            double ini = 0.5;
            // 
            double c = f(ini, u1, un, y, rho);
            if(c > 0){
                // in (0,0.5)
                while(c>0){
                    ini /= 2.0;
                    c = f(ini, u1, un, y, rho);
                }
                return newton(ini, u1, un, y, rho);
            }
            else if(c < 0){
                // in (0.5,1)
                while(c<0){
                    ini = (1.0+ini)/2.0;
                    c = f(ini, u1, un, y, rho);
                }
                return newton(ini, u1, un, y, rho);
            }
            else return ini;
        }

        vector<double> ADMM(vector<double>& Un, int max_iter, double tolerance){
            double rho = 1.0;// penalty function:second order;
            double tau = 2.0;// step_size of Y
            double mu = 5;// dymanic adjustment
            double gamma_p = 2;
            double gamma_d = 2;

            auto U_1 = Un;
            auto U_2 = Un;
            auto U_2_new = Un;

            
            vector<double> Y(N, 0.0);
            
            for(int k=1;k<=max_iter;k++){
                // U1
                auto b = rhs(Un, Y, rho, U_2);
                U_1 = CG(b, U_1, 1e-8, rho, 1e5);
                // U2 Y r s
                double r = 0.0;//primal feasibility
                double s = 0.0;//dual feasibility
                for(int i=0;i<Nx;i++){
                    for(int j=0;j<Ny;j++){
                        int idx = i*Ny+j;
                        // Approximation way
                        // U_2_new[idx] = 1.0-1.0/(1+exp(Y[idx]-theta*(1-2*Un[idx])+rho*(U_1[idx]-U_2_new[idx])));
                        // U_2_new[idx] = (U_2_new[idx]+U_2_old[idx])/2.0;


                        // Newton way
                        U_2_new[idx] = solveu2(U_2[idx], U_1[idx], Un[idx], Y[idx], rho);
                        Y[idx] = Y[idx] + tau*rho*(U_1[idx]-U_2_new[idx]);
                        r += (U_1[idx]-U_2_new[idx])*(U_1[idx]-U_2_new[idx]);
                        s += (U_2[idx]-U_2_new[idx])*(U_2[idx]-U_2_new[idx]);
                    }
                }
                
                U_2 = U_2_new;
                r=sqrt(r);
                s=sqrt(s);
                // cout<<max(r,s)<<endl;
                if(max(r,s) < tolerance) {
                    Un = U_2;
                    cout<<"ADMM converge at "<<k<<endl;
                    break;
                }
                // dymanic adjustment of rho
                /*
                if(r>mu*s){
                    rho = gamma_p*rho;
                }
                else if(s>mu*r){
                    rho = rho/gamma_d;
                }
                else rho = rho;
                rho=max(min(1.0,rho),6.0);

                if(k==max_iter){
                    Un = U_2;
                    cout<<"ADMM not converge."<<endl;
                }
                */
            }
            return Un;
        }
        void solve(){
            for(int n=0;n<Nt;n++){
                auto Un = u[n];
                energy[n]=Energy(Un);
                cout<<energy[n]<<endl;
                u[n+1] = ADMM(Un,1e6,1e-3);
            }
            energy[Nt]=Energy(u[Nt]);
            cout<<energy[Nt]<<endl;
        }
        vector<vector<double>> vectorToarray(const vector<double>& a){
            vector<vector<double>> A(Nx, vector<double>(Ny,0.0));
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    int idx = i*Ny+j;
                    A[i][j] = a[idx];
                }
            }
            return A;
        }
        const vector<vector<double>>& getu() const{
            return u;
        }
        vector<vector<vector<double>>>& getU(){
            for(int n = 0;n<=Nt;n++){
                U[n] = vectorToarray(u[n]);
            }
            return U;
        }
        
};

int main(){
    double dt = 1e10;
    int Nx = 100;
    int Ny = 100;
    int Nt = 30;
    double ep = 0.05;
    AC ACu(dt, Nt, Nx, Ny, ep);
    ACu.solve();
    auto u = ACu.getu();

    auto U =ACu.getU();
    saveDataToFile(U, "u_test_data.txt");
    system("pause");
    return 0;
}