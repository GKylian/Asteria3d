// Asteria3d.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>

#include "Defs/Definitions.h"
#include "Defs/Arrays.h"
#include "prob.h"
#include "Export/vtk.h"
#include "Export/csv.h"

#ifdef PRIMBOUNDS
  #include "Bounds/primBounds.h"
#else
  #include "Bounds/consBounds.h"
#endif // PRIMBOUNDS

#include "Reconstruction/plm.h"
#include "Reconstruction/pcm.h"

#ifdef HLLE
  #include "Riemann/hlle.h"
#endif // HLLE
#ifdef HLLC
  #include "Riemann/hllc.h"
#endif // HLLC
#ifdef HLLD
  #include "Riemann/hlld.h"
#endif // HLLD

#ifdef CT
  #include "Eq/CT.h"
#endif // CT


Arrays u;
Arrays _u;
map2d params;

using namespace std;

bool init(string fname); bool Loop();
ld getdt(ld CFL);
vector<Export> out;


static void  show_usage(string name) {
    cerr << "Usage: " << name << " <option(s)>\n"
        << "Options:\n"
        << "\t-p, --prob\t\tLaunch the simulation with the given file as parameters\n"
        << "\t-h, --help\t\tShow this help message\n"
        << "\t-c, --config\t\tShow the current configuration of Asteria"
        << endl << endl;
}

static void print_setup() {
    cout << "Asteria is configured with the following parameters:" << endl;
#ifdef PNAME
    cout << "\tProblem: \t\t" << PNAME << "\n";
#else
    cout << "\tProblem: \t\tUNKNOWN -> did you specify its name in the problem header (#define PNAME probName) ?\n";
#endif // PNAME

#ifdef MHD
    cout << "\tMHD: \t\t\tyes\n";
#else
    cout << "\tMHD: \t\t\tno\n";
#endif // MHD

#ifdef HLLE
    cout << "\tRiemann solver: \tHLLE\n";
#endif // HLLE
#ifdef HLLC
    cout << "\tRiemann solver: \tHLLC\n";
#endif // HLLC
#ifdef HLLD
    cout << "\tRiemann solver: \tHLLD\n";
#endif // HLLD

#ifdef CHAR
    cout << "\tReconstruction: \tPerform monotonicity constraints in CHARACTERISTIC variables\n";
#else
    cout << "\tReconstruction: \tPerform monotonicity constraints in PRIMITIVE variables\n";
#endif // CHAR

#ifdef TRACING
    cout << "\tChar. tracing: \t\tyes\n";
#else
    cout << "\tChar. tracing: \t\tno\n";
#endif // TRACING




    cout << endl;
}




string pIf(map<string, string> *m, string par) {
    if (m->find(par) != m->end()) {
        return m->at(par);
    }

    return "";
}

bool loadExport() {
    int i = 0;
    for (int i = 0; i < 100; i++)
    {
        string header = "out:"+to_string(i);
        if (params.find(header) != params.end()) {
            cout << "Loading export parameters from header " << header << endl;
            
            map<string, string> p = params[header];
            Export exp;
            if (pIf(&p, "format") == "csv") exp.type = 1;
            else if (pIf(&p, "format") == "vtk") exp.type = 0;
            else { cout << "\tInvalid export format " << pIf(&p, "format") << endl; continue; }

            if (pIf(&p, "dt") != "") exp.dt = stold(pIf(&p, "dt"));
            else { cout << "\tDid not specify the time interval between saves." << endl; continue; }

            for (int v = 0; v < 100; v++) {
                string var = "var"+to_string(v);
                if (pIf(&p, var) != "") {
                    if (strToExp.find(pIf(&p, var)) != strToExp.end()) {
                        exp.exp.push_back(strToExp[pIf(&p, var)]);
                    }
                    else {
                        cout << "\tInvalid variable " << var << " = " << pIf(&p, var) << endl;
                    }
                }
            }
            if (exp.exp.size() == 0) { cout << "\tVector of variables to export is empty !" << endl; continue; }
            
            if (pIf(&p, "ghosts") == "true") exp.ghosts = true;
            else if (pIf(&p, "ghosts") == "false") exp.ghosts = false;
            else if (pIf(&p, "ghosts") == "") exp.ghosts = false;
            else { cout << "\tInvalid value for parameter 'ghosts'. Should be 'true' or 'false'." << endl; continue; }

            if (pIf(&p, "name") != "") exp.name = pIf(&p, "name");
            else { cout << "\tNo name specified, using the id as name" << endl; exp.name = to_string(i); }
            
            out.push_back(exp);
        }


        
    }
    
    cout << endl;
    return true;
}







int main(int argc, char *argv[])
{

    string fname = ""; //The .params file path
    if (argc == 0) {
        show_usage(argv[0]);
        return 1;
    }
    std::vector<std::string> sources; std::string destination;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        }
        else if ((arg == "-p") || (arg == "--prob")) {
            if (i+1 < argc)
                destination = argv[i++];
            else {
                cerr << "ERROR::: -p and --prob require one argument (the .param file path)." << endl << endl;
                return 1;
            }
            fname = argv[i];
        }
        else if ((arg == "-c") || (arg == "--config")) {
            print_setup();
            return 1;
        }
        else {
            show_usage(argv[0]);
            return 1;
        }
    }
    cout << endl;

    if (fname == "") {
        cerr << "ERROR::: No .params file was given. You need to specify a problem file when running Asteria (e.g. ./Asteria -p params/sod.h)." << endl << endl;
        return 1;
    }
    ifstream f(fname);
    if (!f.good()) {
        cerr << "ERROR::: Could not find/open the specified .params file (" << fname << ")." << endl << endl;
        return 1;
    }
    f.close();


    cout << endl << endl;
    /*  Initialize the arrays and set initial conditions  */
    if (!init(fname)) {
        cout << "Asteria3d::main:: The initialization of the simulation failed !" << endl;
        return -1;
    }
    cout << endl << endl;
    

    /* Big loop */
    if (!Loop()) {
        cout << "Asteria3d::main:: The computation failed !" << endl;
    }


    
}



bool init(string fname) {
    /* Coordinates:
    <-->
    ---> Cartesian coordinates: (x,y,z) = (x,y,z)
    ---> Cylindrical coordinates: (x,y,z) = (r,phi,z)
    <-->
    */
    //For cartesian coordinates: (x,y,z) = (x,y,z)

    /* Initialize the arrays with the problem chosen */
    if (!Problem(&u, &out, &params, fname)) {
        cout << "Asteria3D::init:: The initialization of the values from Problem() failed !" << endl; return false;
    }; cout << "Initialized the problem-specific intial values with Problem()." << endl;

    

    /* Read export parameters from 2d map */
    if (!loadExport()) {
        cout << "Asteria3d::init:: Could not load the export parameters from 2d map !" << endl; return false;
    }
    cout << "Loaded the export parameters from 2d map" << endl;

    
    /* Apply the boundary conditions */
    if (!ApplyBounds(&u)) {
        cout << "Asteria3d::init:: Could not apply the boundary conditions !" << endl; return false;
    }
    _u = u;

    /* Temporary: exports the initial data */
    /*exports exp[2] = { exports::PRIM, exports::DIVB };
    vtk_3d(&u, false, exp, 2);*/
    vtk_3d(&u, out);
    csv(&u, &out);

    //DoInLoop(&u); 
    if (!check(&u)) {
        cout << "ERROR:::INIT:: Call to check() after initialization returned false -> aborting the simulation..." << endl;
        return false;
    }

    ld maxC[NVAL] = { -1000 }; ld minC[NVAL] = { 1000 }; fill_n(maxC, NVAL, -1000); fill_n(minC, NVAL, 1000);
    ld maxP[NVAL] = { -1000 }; ld minP[NVAL] = { 1000 }; fill_n(maxP, NVAL, -1000); fill_n(minP, NVAL, 1000);
    
    for(int i = u.i_dl; i <= u.i_dr; i++)
    for(int j = u.j_dl; j <= u.j_dr; j++)
    for(int k = u.k_dl; k <= u.k_dr; k++)
    {
        for (int n = 0; n < NVAL; n++) {

            maxC[n] = fmaxl(maxC[n], u.uC(n, i, j, k));
            minC[n] = fminl(minC[n], u.uC(n, i, j, k));
            maxP[n] = fmaxl(maxP[n], u.uP(n, i, j, k));
            minP[n] = fminl(minP[n], u.uP(n, i, j, k));
        }

    }
    cout << "LOOP:: Maximum values of uC: [" << maxC[0]; for (int n = 1; n < NVAL; n++) cout << ", " << maxC[n]; cout << "]." << endl;
    cout << "LOOP:: Minimum values of uC: [" << minC[0]; for (int n = 1; n < NVAL; n++) cout << ", " << minC[n]; cout << "]." << endl;
    cout << "LOOP:: Maximum values of uP: [" << maxP[0]; for (int n = 1; n < NVAL; n++) cout << ", " << maxP[n]; cout << "]." << endl;
    cout << "LOOP:: Minimum values of uP: [" << minP[0]; for (int n = 1; n < NVAL; n++) cout << ", " << minP[n]; cout << "]." << endl;

    /* Compute the initial dt */
    //cout << test(0.23) << endl;
    ld cfl = stod(params["simulation"]["CFL"]);
    cout << "CFL number: " << cfl << endl;

    ld dt0 = getdt(cfl);
    cout << "INIT:: Initial dt: " << dt0 << endl;
    u.dt = dt0;

    
    
    return true;
}















bool Loop() {
    /*
    ---> To do in loop:
    --->   - Reconstruct values at interfaces from cell averages (don't forget source terms)
    --->   - Compute fluxes at interfaces
    --->   - If MHD: compute EMFs at cell corners and use them to update Bx, By and Bz
    --->   - Update rho to n+1/2
    --->   - Add source terms
    --->   - Update all the variables to n+1
    --->   - Apply boundary conditions
    --->   - Compute a new dt
    --->   - When necessary, export data (.csv, .vtk, ...)
    --->   - Update t
    */
    ld t0 = stod(params["simulation"]["t0"]);
    u.tn = stod(params["simulation"]["tn"]);
    ld cfl = stod(params["simulation"]["CFL"]);
    ld gamma = u.gamma;

    bool nx = u.Nx>1; bool ny = u.Ny>1; bool nz = u.Nz>1; ld hrho = 0.0;

    cout << endl << "LOOP:: Starting computation !" << endl;
    auto start = std::chrono::system_clock::now(); std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    for (u.t = t0; u.t <= u.tn; u.t=u.t+u.dt) {

        //if (u.s == 20) break;
        u.s = u.s+1;
        cout << endl << endl;
        std::cout << "LOOP:: Computing at time " << u.t << " with time-step " << u.dt << endl;
        
        /* ----- 1. Reconstruct interfaces from cell averages ----- */
        cout << "LOOP:: Starting the reconstruction..." << endl;
        for(int i = u.i_dl+1*nx; i <= u.i_dr-1*nx; i++)
        for(int j = u.j_dl+1*ny; j <= u.j_dr-1*ny; j++)
        for(int k = u.k_dl+1*nz; k <= u.k_dr-1*nz; k++)
        {
            
            if (!PLM(&u, i, j, k)) {
                cout << "ERROR:::LOOP:: Reconstruction failed !" << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
        }
        cout << "LOOP:: Reconstructed interfaces from cell averages." << endl;
        

        /* ----- 2. Compute fluxes from newly computed interface values ----- */
        cout << "LOOP:: Starting to compute the fluxes..." << endl;

        
        for(int i = u.i_dl+2*nx; i <= u.i_dr-1*nx; i++)
        for(int j = u.j_dl+2*ny; j <= u.j_dr-1*ny; j++)
        for(int k = u.k_dl+2*nz; k <= u.k_dr-1*nz; k++)
        {
            if (!getFluxes(&u, i, j, k)) {
                cout << "ERROR:::LOOP:: getFluxes failed !" << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
            //TODO:cout << "\t->" << u.Fy(NVAL-2, i, j, k) << endl;
            for (int n = 0; n < NVAL; n++) {
                if (isnan(u.Fx(n, i, j, k))) {
                    cout << "ERROR:::LOOP:: Fx[" << n << "] is NaN !" << endl; return false;
                }
                if (isnan(u.Fy(n, i, j, k))) {
                    cout << "ERROR:::LOOP:: Fy[" << n << "] is NaN !" << endl; return false;
                }
                if (isnan(u.Fz(n, i, j, k))) {
                    cout << "ERROR:::LOOP:: Fz[" << n << "] is NaN !" << endl; return false;
                }
                
                    
            }

        }
        cout << "LOOP:: Computed fluxes." << endl;
        

        /* ----- 4. Update to n+1 ----- */
        for(int i = u.i_dl+3*nx; i <= u.i_dr-2*nx; i++)
        for(int j = u.j_dl+3*ny; j <= u.j_dr-2*ny; j++)
        for(int k = u.k_dl+3*nz; k <= u.k_dr-2*nz; k++)
        {
            /* --- 4.1) Compute rho at t(n+1/2) --- */
            hrho = u.uC(0, i, j, k);
            if (u.Nx > 1) hrho += -u.dt/(2.0*u.dx) * (u.Fx(0, i+1, j, k)-u.Fx(0, i, j, k));
            if (u.Ny > 1) hrho += -u.dt/(2.0*u.dy) * (u.Fy(0, i, j+1, k)-u.Fy(0, i, j, k));
            if (u.Nz > 1) hrho += -u.dt/(2.0*u.dz) * (u.Fz(0, i, j, k+1)-u.Fz(0, i, j, k));


            /* --- 4.2) Add the gravitational source terms to _u --- */
            
            /* Gravity with acceleration \vec{g} */
            if (gx != 0.0 || gy != 0.0 || gz != 0.0) {      //TODO: '-=' or '+=' ?
                _u.uC(1, i, j, k) -= u.dt*hrho*gx;
                _u.uC(2, i, j, k) -= u.dt*hrho*gy; 
                _u.uC(3, i, j, k) -= u.dt*hrho*gz;
#ifndef ISO
                if (u.Nx > 1) _u.uC(4, i, j, k) -= u.dt*gx*(u.Fx(0, i+1, j, k)+u.Fx(0, i, j, k))/2.0;
                if (u.Ny > 1) _u.uC(4, i, j, k) -= u.dt*gy*(u.Fy(0, i, j+1, k)+u.Fy(0, i, j, k))/2.0;
                if (u.Nz > 1) _u.uC(4, i, j, k) -= u.dt*gz*(u.Fz(0, i, j, k+1)+u.Fz(0, i, j, k))/2.0;
#endif // !ISO

            }
            /* Gravity with static potential defined in prob.h */
            else if (Phii(&u, i, j, k) != 0.0){
                ld x[3] = { 0 }; u.pos(x, i, j, k); ld Phic = Phii(&u, i, j, k); ld Phir, Phil = 0.0;
                if (u.Nx>1) {
                    Phir = Phi(x[0]+0.5*u.dx, x[1], x[2]); Phil = Phi(x[0]-0.5*u.dx, x[1], x[2]);
                    _u.uC(1, i, j, k) -= u.dt/u.dx*hrho*(Phir-Phil);
#ifndef ISO
                    _u.uC(4, i, j, k) -= u.dt/u.dx*(u.Fx(0, i, j, k)*(Phic-Phil) + u.Fx(0, i+1, j, k)*(Phir-Phic));
#endif // !ISO
                }
                if (u.Ny>1) {
                    Phir = Phi(x[0], x[1]+0.5*u.dy, x[2]); Phil = Phi(x[0], x[1]-0.5*u.dy, x[2]);
                    _u.uC(2, i, j, k) -= u.dt/u.dy*hrho*(Phir-Phil);
#ifndef ISO
                    _u.uC(4, i, j, k) -= u.dt/u.dy*(u.Fy(0, i, j, k)*(Phic-Phil) + u.Fy(0, i, j+1, k)*(Phir-Phic));
#endif // !ISO
                }
                if (u.Nz>1) {
                    Phir = Phi(x[0], x[1], x[2]+0.5*u.dz); Phil = Phi(x[0], x[1], x[2]-0.5*u.dz);
                    _u.uC(3, i, j, k) -= u.dt/u.dz*hrho*(Phir-Phil);
#ifndef ISO
                    _u.uC(4, i, j, k) -= u.dt/u.dz*(u.Fz(0, i, j, k)*(Phic-Phil) + u.Fz(0, i, j, k+1)*(Phir-Phic));
#endif // !ISO
                }

            }




            /* --- 4.3) Update u to n+1 in _u using the fluxes previously computed --- */
            bool BxCT = false, ByCT = false, BzCT = false; /* Have Bx, By and/or Bz been updated with CT ? */
#ifdef CT
            ld dE = 0.0; /* Stores the differences every time to check if very small */
            if (u.Ny>1 && u.Nz>1) { /* Update By and Bz using Ey */
                dE = -u.dt/u.dz * (getEx(&u, i, j, k+1)-getEx(&u, i, j, k)); if (null(dE)) dE = 0.0;
                _u.uC(NVAL-2, i, j, k) += dE; ByCT = true;
                dE = u.dt/u.dy * (getEx(&u, i, j+1, k)-getEx(&u, i, j, k)); if (null(dE)) dE = 0.0;
                _u.uC(NVAL-1, i, j, k) += dE; BzCT = true;
                
            }
            if (u.Nx>1 && u.Nz>1) { /* Update Bx and Bz using Ey */
                dE = u.dt/u.dz * (getEy(&u, i, j, k+1)-getEy(&u, i, j, k)); if (null(dE)) dE = 0.0;
                _u.uC(NVAL-3, i, j, k) += dE; BxCT = true;
                dE = -u.dt/u.dx * (getEy(&u, i+1, j, k)-getEy(&u, i, j, k)); if (null(dE)) dE = 0.0;
                _u.uC(NVAL-1, i, j, k) += dE; BzCT = true;
            }
            if (u.Nx > 1 && u.Ny > 1) { /* Update Bx and By using Ez */
                dE = -u.dt/u.dy * (getEz(&u, i, j+1, k)-getEz(&u, i, j, k)); if (null(dE)) dE = 0.0;
                _u.uC(NVAL-3, i, j, k) += dE; BxCT = true;
                dE = u.dt/u.dx * (getEz(&u, i+1, j, k)-getEz(&u, i, j, k)); if (null(dE)) dE = 0.0;
                _u.uC(NVAL-2, i, j, k) += dE; ByCT = true;
            }
            
#endif // CT

            if (u.Nx>1 && u.dim==1 && (u.Fx(NVAL-3, i+1, j, k)-u.Fx(NVAL-3, i, j, k)) != 0.0) {
                cout << "\tThe domain is (N,0,0), but dFx[Bx] is not null !" << endl;
                cout << "\t Fx[Bx](i+1) = " << u.Fx(NVAL-3, i+1, j, k) << ", Fx[Bx](i) = " << u.Fx(NVAL-3, i, j, k) << endl;
                return false;
            }
            if (u.Ny>1 && u.dim==1 && (u.Fy(NVAL-2, i, j+1, k)-u.Fy(NVAL-2, i, j, k)) != 0.0) {
                cout << "\tThe domain is (0,N,0), but dFy[By] is not null !" << endl;
                cout << "\t\tAt (" << i << ", " << j << ", " << k << ")." << endl;
                cout << "\t\tj+1 = " << u.Fy(NVAL-2, i, j+1, k) << ",   j = " << u.Fy(NVAL-2, i, j, k) << endl;
                return false;
            }

            for (int n = 0; n < NVAL; n++) {
#ifdef MHD
                if (n == NVAL-3 && BxCT == true) continue; /* Skip if the component of B has already */
                if (n == NVAL-2 && ByCT == true) continue; /* been updated using the CT algorithm    */
                if (n == NVAL-1 && BzCT == true) continue;
#endif // MHD
                
                ld dFxdx = 0.0; if (u.Nx>1) dFxdx = u.dt/u.dx * (u.Fx(n, i+1, j, k)-u.Fx(n, i, j, k));
                ld dFydy = 0.0; if (u.Ny>1) dFydy = u.dt/u.dy * (u.Fy(n, i, j+1, k)-u.Fy(n, i, j, k));
                ld dFzdz = 0.0; if (u.Nz>1) dFzdz = u.dt/u.dz * (u.Fz(n, i, j, k+1)-u.Fz(n, i, j, k));
                _u.uC(n, i, j, k) += -dFxdx - dFydy - dFzdz;

                //if (n == 7 /*&& _u.uC(n, i, j, k) == 0.0*/) cout << "Bz == 0.0 !" << endl;

                /*if (null(dFxdx) && null(dFydy) && null(dFydy)) {
                    cout << "_u.uC[" << n << "] doesn't change !" << endl;
                    cout << "---> dFxdx: " << dFxdx << ", dFydy: " << dFydy << ", dFzdz: " << dFzdz << endl;
                    cout << "---> F(i+1) = " << u.Fx(n, i+1, j, k) << ", f(i) = " << u.Fx(n, i, j, k) << endl;
                }*/

                if (isnan(_u.uC(n, i, j, k))) {
                    cout << "ERROR:::LOOP:: After evolution to n+1, uC[" << n << "] is NaN !" << endl;
                    cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                    cout << "---> dFndn = " << dFxdx << ", " << dFydy << ", " << dFzdz << endl;
                    cout << "---> dFn = " << (u.Fx(n, i+1, j, k)-u.Fx(n, i, j, k)) << ", " << (u.Fy(n, i, j+1, k)-u.Fy(n, i, j, k)) << ", " << (u.Fz(n, i, j, k+1)-u.Fz(n, i, j, k)) << endl;
                    cout << "---> Fn(n+1) = " << u.Fx(n, i+1, j, k) << ", " << u.Fy(n, i, j+1, k) << ", " << u.Fz(n, i, j, k+1) << endl;
                    return false;
                }
            }
        }
        for(int i = u.i_dl; i <= u.i_dr; i++)
        for(int j = u.j_dl; j <= u.j_dr; j++)
        for(int k = u.k_dl; k <= u.k_dr; k++)
        {
            for (int n = 0; n < NVAL; n++) {
                u.uC(n, i, j, k) = _u.uC(n, i, j, k);
                //u.uP(n, i, j, k) = _u.uP(n, i, j, k);
            }
        }
        cout << "LOOP:: Updated to n+1." << endl;
        

        /* ----- 5. Get primitive variables and B center ----- */
        for(int i = u.i_dl; i <= u.i_dr; i++)
        for(int j = u.j_dl; j <= u.j_dr; j++)
        for(int k = u.k_dl; k <= u.k_dr; k++)
        {
#ifdef MHD
            BCenter(&u, i, j, k);
#endif

            toPrimitive(&u, i, j, k, gamma);

            
        }



        /* ----- 6. Apply boundary conditions ----- */
        ApplyBounds(&u); //TODO: Transfer the handling of the BCs to the problem header.


        /* ----- 7. Call DoInLoop in the problem header ----- */
        if(!DoInLoop(&u)) return false;


        /* ----- 8. Re-apply boundary conditions in case we changed values in DoInLoop ----- */
        ApplyBounds(&u); //TODO: Transfer the handling of the BCs to the problem header.
        //Also check max/min values of uP and uC
        ld maxC[NVAL] = { -1000 }; ld minC[NVAL] = { 1000 };
        ld maxP[NVAL] = { -1000 }; ld minP[NVAL] = { 1000 };
        for(int i = u.i_dl; i <= u.i_dr; i++)
        for(int j = u.j_dl; j <= u.j_dr; j++)
        for(int k = u.k_dl; k <= u.k_dr; k++)
        {
            for (int n = 0; n < NVAL; n++) {
                _u.uC(n, i, j, k) = u.uC(n, i, j, k);
                _u.uP(n, i, j, k) = u.uP(n, i, j, k);

                maxC[n] = fmaxl(maxC[n], u.uC(n, i, j, k));
                minC[n] = fminl(minC[n], u.uC(n, i, j, k));
                maxP[n] = fmaxl(maxP[n], u.uP(n, i, j, k));
                minP[n] = fminl(minP[n], u.uP(n, i, j, k));


                if (u.uC(0, i, j, k) <= 0.0 || null(u.uC(0, i, j, k))) {
                cout << "ERROR:::LOOP:: rho in uC is either NaN or null/negative: " << u.uC(0, i, j, k) << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
            if (u.uP(0, i, j, k) <= 0.0 || null(u.uP(0, i, j, k))) {
                cout << "ERROR:::LOOP:: rho in uP is either NaN or null/negative: " << u.uP(0, i, j, k) << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }

            if (u.uC(4, i, j, k) <= 0.0 || null(u.uC(4, i, j, k))) {
                cout << "ERROR:::LOOP:: E in uC is either NaN or null/negative: " << u.uC(4, i, j, k) << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
            if (u.uP(4, i, j, k) <= 0.0 || null(u.uP(4, i, j, k))) {
                cout << "ERROR:::LOOP:: P in uP is either NaN or null/negative: " << u.uP(4, i, j, k) << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
            }
        }
        cout << "LOOP:: Maximum values of uC: "; pArray(maxC, NVAL); cout << "." << endl;
        cout << "LOOP:: Minimum values of uC: "; pArray(minC, NVAL); cout << "." << endl;
        cout << "LOOP:: Maximum values of uP: "; pArray(maxP, NVAL); cout << "." << endl;
        cout << "LOOP:: Minimum values of uP: "; pArray(minP, NVAL); cout << "." << endl;

        /* ----- 9. Call the Check() function in the problem file ----- */
        if (!check(&u)) {
            cout << "ERROR:::LOOP:: Call to check() returned false -> stoping the loop..." << endl;
            return false;
        }

        /* ----- 10. Compute a new dt ----- */
        u.dt = getdt(cfl);


        
        /* ----- 11. Export (if necessary) the variables ----- */
        vtk_3d(&u, out);
        csv(&u, &out);


    }
    vtk_3d(&u, out);
    csv(&u, &out);

    auto end = std::chrono::system_clock::now(); std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::chrono::duration<double> elapsed = end-start;
    cout << endl << "LOOP:: Finished the computation !" << endl;

    cout << "\tStarted at time " << std::ctime(&start_time);
    cout << "\tFinished at time " << std::ctime(&end_time);
    cout << "\tElapsed time (in s): " << elapsed.count() << endl;

    return true;
}



ld getdt(ld CFL) {
    ld dt = 10000.0;
    
    ld cx[WV] = { 0 }; ld cy[WV] = { 0 }; ld cz[WV] = { 0 };
    ld dts[3] = { 0 };
    ld Cx = 0.0, Cy = 0.0, Cz = 0.0;
    ld gamma = u.gamma;

    // Loop through the whole domain (including ghost cells)
    for(int i = u.i_dl; i <= u.i_dr; i++)
    for(int j = u.j_dl; j <= u.j_dr; j++)
    for(int k = u.k_dl; k <= u.k_dr; k++)
    {
        
        getWavespeeds(&u, i, j, k, gamma, cx, cy, cz); 
        Cx = (ld)*std::max_element(cx, cx+WV); Cy = (ld)*std::max_element(cy, cy+WV); Cz = (ld)*std::max_element(cz, cz+WV);
        

        dts[0] = u.dx/(u.uP(1, i, j, k)+Cx); if (u.Nx<=1) dts[0] = 10000.0; if (isnan(dts[0])) dts[0] = 10000.0;
        dts[1] = u.dy/(u.uP(2, i, j, k)+Cy); if (u.Ny<=1) dts[1] = 10000.0; if (isnan(dts[1])) dts[1] = 10000.0;
        dts[2] = u.dz/(u.uP(3, i, j, k)+Cz); if (u.Nz<=1) dts[2] = 10000.0; if (isnan(dts[2])) dts[2] = 10000.0;
        ld thdt = CFL * (*std::min_element(dts, dts+3));

        if (thdt <= 0) {
            thdt = 100.0;
        }

        dt = fminl(dt, thdt);

    }
    cout << "Computed dt: " << dt << endl;

    return dt;

}






