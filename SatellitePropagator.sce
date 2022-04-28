function Initialise(Satellite,kep,cjd,wi,Euleri)
/*----------------------------------------------------------------------
this function helps load the geometry of the and the initial position 
and orientation of the satellite into the simulation

function parameters:
Satellite: An object describing the satellite geometry and hardware
kep: intial keplerian orbital elements
cjd: initial modified julian time
wi: initial angular velocity
Euleri: initial Euler Angles
----------------------------------------------------------------------*/    
    
    //Unpacks the geometry of the satellite
    A = Satellite.A'
    B = Satellite.B'
    C = Satellite.C'
   
    //inital position and velocity of satellite
    [Ri, vi] = CL_oe_kep2car(kep)
    
    //transform from body frame to initial orbital frame
    RMatrixi = CL_fr_lvlhMat(Ri, vi)'
    
    /*Shifts the position of the satellite's COM to the origin so that
    rotations can be more easily applied*/                 
    COM = Satellite.COM' 
    N = size(A)(1,2)
    COM = repmat(COM,1,N)
    A = A - COM
    B = B - COM
    C = C - COM
    
    /*checks if satellite has solarpanels, initializes solarpanels at
    origin and calculates the direction vector of the rotation axis*/
    if typeof(Satellite.solarpanels) == "st" then
        D = Satellite.solarpanels.D'
        E = Satellite.solarpanels.E'
        Nsolarpanel = size(D)(1,2)
        COM = Satellite.COM' 
        COM = repmat(COM,1,Nsolarpanel)
        D = D - COM
        E = E - COM
        k = (E-D)./repmat((sum((E-D).^2,1).^0.5),3,1)
        Satellite.solarpanels.D = D
        Satellite.solarpanels.k = k
    end
     
    //calculates initial geometric parameters
    AB = B - A
    AC = C - A
    OM = A + 0.5*AB + 0.5*AC
    n = cross(AB,AC)
    An = (sum(AB.^2,1).^0.5).*(sum(AC.^2,1).^0.5)
       
    //initial atmospheric conditions
    [dens,Ti,M] = NRLMSISE00(cjd,CL_co_car2ell(Ri))
        
    //convert vectors to unit vector
    AB = AB./repmat(sum(AB.^2,1).^0.5,3,1)
    AC = AC./repmat(sum(AC.^2,1).^0.5,3,1)    
    n = n./repmat(sum(n.^2,1).^0.5,3,1)
    Aref = sum(abs([1,0,0]*n.*An),2)/2
    
    //packaging parameters
    Satellite.A = A
    Satellite.B = B
    Satellite.C = C
    Satellite.OM = OM
    Satellite.AB = AB
    Satellite.AC = AC
    Satellite.An = An
    Satellite.n = n
    Satellite.Aref = Aref
    Satellite.initial_coordinates = Ri
    Satellite.initial_velocity = vi
    Satellite.initial_orientation = RMatrixi
    Satellite.dens = dens
    Satellite.Ti = Ti
    Satellite.M = M       
    Ri.values = Ri'
    Ri.time = 0
    vi.values = vi'
    vi.time = 0
    cjd.values = cjd
    cjd.time = 0
    wi.values = wi'
    wi.time = 0
    qi = CL_rot_angles2quat([3,2,1], Euleri)
    qi = [real(qi);imag(qi)]
    qi.values = qi'/norm(qi)
    qi.time = 0
    [Ri,vi,cjd,wi,qi,kep,Satellite] = return(Ri,vi,cjd,wi,qi,kep,Satellite)
endfunction

function Satellite = GeometryUpdater(Satellite,q,R,v)
/*----------------------------------------------------------------------
this function updates the geometry of the satellite by applying the 
rotation matrix onto the geomtric parameters of the satellite
----------------------------------------------------------------------*/ 
    
    /*finds the rotation matrix describing the transformation between
    the body frame of the satellite and the ECI frame*/
    q = CL_rot_defQuat(q)  
    RMatrix = Satellite.initial_orientation*CL_rot_quat2matrix(q)'
    
    //rotates body of satellite
    Satellite.OM = RMatrix*Satellite.OM
    Satellite.AB = RMatrix*Satellite.AB
    Satellite.AC = RMatrix*Satellite.AC
    Satellite.n = RMatrix*Satellite.n
    
    //rotates the solarpanels if there are any
    if typeof(Satellite.solarpanels) == "st" then
        Satellite.solarpanels.D = RMatrix*Satellite.solarpanels.D
        Satellite.solarpanels.k = RMatrix*Satellite.solarpanels.k
    end

    Satellite.RMatrix = RMatrix
endfunction

function [Satellite,Solarflux,SolarPanelAngle] = solarpower(Satellite,pos_sun,R)
/*----------------------------------------------------------------------
this function helps to futher rotate the solar panels such that the
panels recieve the maximum amount of solar flux
----------------------------------------------------------------------*/

    //checks for the presence of solar panels     
    if typeof(Satellite.solarpanels) == "boolean" & Satellite.solarpanels == %f then
        [Satellite,Solarflux,SolarPanelAngle] = return(Satellite,0,0)        
    else
        distance = norm(pos_sun - R) //distance from sun to satellite
        alignmentvector = (pos_sun - R) //unit direction vector to sun
        
        //geometric paramters of the solar panels
        idxtop = Satellite.solarpanels.idxtop
        idxbot = Satellite.solarpanels.idxbot
        n = Satellite.n(:,idxtop)
        An = Satellite.An(idxtop)
        
        //check if the solarpanels are set to be aligned to the sun
        if Satellite.solarpanels.align == %t then        
            OM = Satellite.OM(:,idxtop)
            AB = Satellite.AB(:,idxtop)
            AC = Satellite.AC(:,idxtop)        
            D = Satellite.solarpanels.D
            k = Satellite.solarpanels.k
            
            //calculates optimal angle of rotation
            x = cross(k,n)
            x = alignmentvector'*x/norm(x)
            z = alignmentvector'*n
            t = atan(x./z)
            
            //rotates the panels through quaternions 
            spanelN =  size(idxtop)(1,2)    
            OM = OM - D
            qOM =  CL_rot_defQuat(zeros(1,spanelN),OM)
            qAB =  CL_rot_defQuat(zeros(1,spanelN),AB) 
            qAC =  CL_rot_defQuat(zeros(1,spanelN),AC) 
            qn =  CL_rot_defQuat(zeros(1,spanelN),n)
            q2 = CL_rot_defQuat(cos(t/2),k.*repmat(sin(t/2),3,1))
            OM = imag(q2*qOM*q2')+D
            AB = imag(q2*qAB*q2')
            AC = imag(q2*qAC*q2')
            n = imag(q2*qn*q2')                     
            Satellite.OM(:,idxtop) = OM
            Satellite.OM(:,idxbot) = OM
            Satellite.n(:,idxtop) = n
            Satellite.n(:,idxbot) = -n
            Satellite.AB(:,idxtop) = AB
            Satellite.AC(:,idxtop) = AC
            Satellite.AB(:,idxbot) = AC
            Satellite.AC(:,idxbot) = AB           
            SolarPanelAngle = t(1,1)
        else 
            SolarPanelAngle = 0    
        end
        
        //calculates incident solarflux on solarpanels
        eclrat = CL_gm_eclipseCheck(R, pos_sun)
        Aref = abs(alignmentvector'*n.*An/norm(n))
        Solarflux = sum(Aref*eclrat*3.864*10^26/(4*%pi*distance^2),2)
    end
endfunction

function [Aero_Acceleration,Aero_Torque,Cd] = PerPlaneForce(Satellite,v,w)
/*----------------------------------------------------------------------
this function is the implementation of sentman's equations for free 
molecular flow aerodynamics. The satellite geometry is passed into the
function to calculate the aerodynamic acceleration and torque
----------------------------------------------------------------------*/
    v = v'
    vcap = v/norm(v)
    Ti = Satellite.Ti
    dens = Satellite.dens
    M = Satellite.M
    
    //temperature of reflected molecules
    Tr = Satellite.surfacetemperature
    
    //unpacks the needed geometric parameters     
    OM = Satellite.OM
    AB = Satellite.AB
    AC = Satellite.AC
    n = Satellite.n
    An = Satellite.An
           
    //find mass velocity components in directions n, AB and AC
    N = size(OM)(1,2)
    vw = cross(repmat(w,1,N),OM) + repmat(v',1,N)
    vn = sum(vw.*n,1)
    vAB = sum(vw.*AB,1)
    vAC = sum(vw.*AC,1)
        
    B = M/(2*8.31*Ti)
    Br = M/(2*8.31*Tr)
    
    //average velocity of particles colliding with the plane in the direction of normal
    Vave = (vn/2).*(1 + erf(vn*sqrt(B))) + (1/(2*sqrt(%pi*B)))*%e^(-B*(vn.**2))
        
    //density of particles leaving the plane    
    dens_r = dens*sqrt(Br)*(vn*sqrt(%pi).*(1 + erf(vn*sqrt(B))) + (1/sqrt(B))*%e^(-B*vn.**2))
        
    //average mean square velocity of particles colliding with the plane
    Vmeansq = (vn/(2*sqrt(%pi))).*(vn*sqrt(%pi).*(1+erf(vn*sqrt(B))) + (1/sqrt(B))*%e^(-B*vn.**2)) + (1/(4*B))*(1+erf(vn*sqrt(B)))
        
    //calculates force in the n, AB and AC directions, 8.31*Tr/M is mean square of reflected particles
    Fn = An.*(dens*Vmeansq + dens_r/(4*Br))
    FAB = dens*vAB.*Vave.*An
    FAC = dens*vAC.*Vave.*An 
        
    //converts Force back to body frame of the satellite
    Fn = repmat(Fn,3,1)
    FAB = repmat(FAB,3,1)
    FAC = repmat(FAC,3,1)
    ForceOnPlane = -(Fn.*n + FAB.*AB + FAC.*AC)
    Force = sum(ForceOnPlane,2)
    Aero_Acceleration = Force/Satellite.Mass

    //Torque calculation
    Torque = cross(OM,ForceOnPlane)
    Aero_Torque = inv(Satellite.RMatrix)*sum(Torque,2)
    
    //Drag coefficient calculation
    Cd = 2*abs(vcap*Force/(dens*Satellite.Aref*norm(v)^2))
endfunction

function [gEarth,gTorque,gSun] = Gravity(R,Satellite,w,pos_sun)
/*----------------------------------------------------------------------
Calculates the gravitational force experienced by the satellite by
spherical harmonics, the third body effect from the sun as well as
the gravity gradient torque due to the Earth's gravity  
----------------------------------------------------------------------*/
    gEarth = CL_fo_sphHarmAcc(R,[3,0],%CL_eqRad,%CL_mu,%CL_j1jn,%CL_cs1nm,%t)
    mu_sun = 1.327*10^20
    gSun = CL_fo_thirdBodyAcc(R, pos_sun, mu_sun)
    gAcceleration = gEarth + gSun
    gTorque = (3*%CL_mu/norm(R)**3)*cross(Satellite.RMatrix*-R,Satellite.MOI*Satellite.RMatrix*-R)
endfunction

function [MTorque,Dipole] = MagneticTorque(B,Satellite,w)
/*----------------------------------------------------------------------
Calculates the Magnetic torque and the resulting dipole moment of three
axis magnetorquers. The algorithm used is the "B-dot" algorithm  
----------------------------------------------------------------------*/
    if typeof(Satellite.Magnetorquer) == "boolean" & Satellite.Magnetorquer == %f then
        [MTorque,Dipole] = return([0;0;0],[0;0;0])
    else
        B = inv(Satellite.RMatrix)*B
        BDot = cross(B,w)
        Dipole = -Satellite.Magnetorquer.gain*BDot
        MTorque = cross(Dipole,B)
    end   
endfunction

function [Acceleration,AngularAccel,Output] = propagator(R,v,q,w,dens,Ti,M,B,cjd)
/*----------------------------------------------------------------------
The superfunction combining many of the previous functions to calculate
the total acceleration and angular acceleration acting on the satellite
to propagate the satellite one time step forward in the simulation  
----------------------------------------------------------------------*/    
    Satellite = Satellite
    Satellite.dens = dens
    Satellite.Ti = Ti
    Satellite.M = M
    pos_sun = CL_eph_sun(cjd)
    q = q/norm(q)
    Satellite = GeometryUpdater(Satellite,q,R,v)   
    [Satellite,Solarflux,SolarPanelAngle] = solarpower(Satellite,pos_sun,R)
    [Aero_Acceleration,Aero_Torque,Cd] = PerPlaneForce(Satellite,v,w)
    [gEarth,gTorque,gSun] = Gravity(R,Satellite,w,pos_sun)
    [MTorque,Dipole] = MagneticTorque(B,Satellite,w)
    Acceleration = Aero_Acceleration + gEarth + gSun
    Torque = Aero_Torque + gTorque + MTorque
    AngularAccel = inv(Satellite.MOI)*(Torque - cross(w,Satellite.MOI*w))
    Output = [R;v;q;w;Dipole;Solarflux;SolarPanelAngle;Cd]
endfunction

function qdot = quatRate(q,w)
/*----------------------------------------------------------------------
Evaluates the quaternion time derivative from the body angular velocity.  
----------------------------------------------------------------------*/    
    q = CL_rot_defQuat(q)
    w = CL_rot_defQuat(0,w)
    qdot = 0.5*q*w
    qdot = [real(qdot);imag(qdot)]
endfunction

function [dens,Ti,M] = NRLMSISE00(cjd,pos_ell)
/*----------------------------------------------------------------------
Implements the NRL-MSISE-00 atmospheric model from celestlab and sets
default f107 and f107a as 150 and ap = 15 
----------------------------------------------------------------------*/
    lon = pos_ell(1,1)
    lat = pos_ell(2,1)
    alt = pos_ell(3,1)
    f107 = 150;
    f107a = 150;
    ap = 15;
    Mr = [4,16,28,32,40,1,14,16]
    [dens, Ti, temp_exo, dens_part] = CL_mod_atmMSIS00(cjd, lon, lat, alt, f107, f107a, ap, res=["dens", "temp", "temp_exo", "dens_part"])
    M = Mr*dens_part/(dens*1000) 
endfunction

function [pos_ell] = ECItoECEF(pos_GCRS,cjd)
/*----------------------------------------------------------------------
Implements the frame conversions from celestlab to convert the cartesian
coordinates in GCRS frame to the elliptical coordinates in the ITRS frame  
----------------------------------------------------------------------*/
    pos_ITRS = CL_fr_convert("GCRS", "ITRS", cjd, pos_GCRS)
    pos_ell = CL_co_car2ell(pos_ITRS)
endfunction



/*----------------------------------------------------------------------
The functions below serve to processs data from the simulation and
present the data visually such as through 2D ground track plots of the 
satellite or through animating the Satellite motion into a gif  
----------------------------------------------------------------------*/

function GroundTrackPlot(lon,lat,colour)
/*----------------------------------------------------------------------
Plots the animated ground track of the satellite

Parameters:
lon: Nx1 matrix of longitude values
lat: Nx1 matrix of the corresponding lattitude values
colour: color id of the plot 
----------------------------------------------------------------------*/
    plot(0,0)
    ax = gca()
    ax.tight_limits=["on","on"]
    ax.data_bounds=[-180 -90;180 90]
    j = [];
    for i = 2:length(lon)
        if ((lon(i) > 170) && (lon(i-1) < -170)) || ((lon(i) < -170) && (lon(i-1) > 170))
            j = [j,i]
        end
    end 
    j = [j,length(lon)]
    comet(lon(1:(j(1)-1)),lat(1:(j(1)-1)),"colors",colour)
    for i = 1:(length(j)-1)
        comet(lon(j(i):(j(i+1)-1)),lat(j(i):(j(i+1)-1)),"colors",colour)
    end 
    xlabel("Longitude (deg)","fontsize",4)
    ylabel("Latitude (deg)","fontsize",4)
endfunction

function OrbitView(R,colour)
/*----------------------------------------------------------------------
Plots the animated orbit of the satellite

Parameters:
R: Cartesian coordinates of the satellite in the ECI frame
colour: color id of the plot  
----------------------------------------------------------------------*/
    xbound = max([abs(R(:,1));6378.e3])
    ybound = max([abs(R(:,2));6378.e3])
    zbound = max([abs(R(:,3));6378.e3])
    plot3d([0,0],[0,0],[0,0])
    ax = gca()
    ax.data_bounds=[-xbound, -ybound, -zbound;xbound, ybound, zbound]
    comet3d(R(:,1),R(:,2),R(:,3),"colors",colour)
endfunction

function DataProcess(Mode,filename,Output,pos_ell)
/*----------------------------------------------------------------------
Processes data from the simulation by saving the data onto a .csv file 
or load data from a .csv file

Paramters:
Mode: "load" to load data from a .csv file 
      "extract" to extract data to a .csv file
      "extractandload" to do both
filename: name of the file to be saved
Output: simulation output 
pos_ell: simulation output of the elliptical coordinates of the satellite

exact filepath can be configured at the start of the function below
----------------------------------------------------------------------*/
    //exact filepath can be configured here
    filepath = "C:\Users\TTzeYoun\Documents\Simulations\"
    select Mode
    case "load" then
        Data = csvRead(filepath+filename+".csv")
        time = Data(2:$,1)
        lon = Data(2:$,2)*180/%pi - 180
        lat = Data(2:$,3)*180/%pi
        alt = Data(2:$,4)
        R = Data(2:$,5:7)
        v = Data(2:$,8:10)
        q =  CL_rot_defQuat(Data(2:$,11:14)')
        w = Data(2:$,15:17)
        Dipole = Data(2:$,18:20)
        Solarflux = Data(2:$,21)
        SolarPanelAngle = Data(2:$,22)
        Cd = Data(2:$,23)
        Ri = Data(1,2:4)
        vi = Data(1,5:7)
        cjd = Data(1,8)
        qi = CL_rot_defQuat(Data(1,9:12)')
        wi = Data(1,13:15)
        kep = Data(1,16:21)' 
        M1 = CL_rot_quat2matrix(q)'
        RMatrixi = CL_fr_lvlhMat(R', v')
        M2 = []
        for i = 1 : size(time)(1,1)
            M2(:,:,i) = RMatrixi(:,:,i)*M1(:,:,i)
        end             
        Angles = CL_rot_matrix2angles(M2,[3,2,1])     
        [time,lon,lat,alt,R,v,q,w,Dipole,Solarflux,SolarPanelAngle,Cd,Ri,vi,cjd,qi,wi,kep,Angles] = return(time,lon,lat,alt,R,v,q,w,Dipole,Solarflux,SolarPanelAngle,Cd,Ri,vi,cjd,qi,wi,kep,Angles)
    case "extract" then
        comments = strcat([Satellite.Name,string([Ri.values,vi.values,cjd.values,qi.values,wi.values,kep',zeros(1,2)])],",")
        csvWrite([Output.time,pos_ell.values,Output.values],filepath+filename+".csv",[],[],[],comments)
    case "extractandload"
        DataProcess("extract",filename,Output,pos_ell)
        DataProcess("load",filename)
        [time,lon,lat,alt,R,v,q,w,Dipole,Solarflux,SolarPanelAngle,Cd,Ri,vi,cjd,qi,wi,kep,Angles] = return(time,lon,lat,alt,R,v,q,w,Dipole,Solarflux,SolarPanelAngle,Cd,Ri,vi,cjd,qi,wi,kep,Angles)               
    end
endfunction

function Animation(Satellite,filename)
/*----------------------------------------------------------------------
Animates the satellite motion from data in the .csv file and saves the
animation into a gif

exact filepath can be configured below 
----------------------------------------------------------------------*/
    DataProcess("load",filename)
    Initialise(Satellite,kep,cjd,wi,CL_rot_quat2angles(qi, [3,2,1]))
    Ai = Satellite.A
    Bi = Satellite.B
    Ci = Satellite.C
    ABi = Bi - Ai
    ACi = Ci - Ai
    Fi = Ai + ABi + ACi 
    N = size(time)(1,1)
    Nplate = size(Ai)(1,2) 
    M = CL_rot_quat2matrix(q)'
    //exact filepath can be configured here
    idGif = animaGIF(gcf(),'C:\Users\TTzeYoun\Documents\Simulations\'+filename+'.gif',20)
    for j = 1 : N      
        RMatrix = Satellite.initial_orientation*M(:,:,j)
        A = RMatrix*Ai
        B = RMatrix*Bi
        C = RMatrix*Ci
        F = RMatrix*Fi
        if typeof(Satellite.solarpanels) == "st" & Satellite.solarpanels.align == %t then
            idxtop = Satellite.solarpanels.idxtop
            idxbot = Satellite.solarpanels.idxbot
            D = RMatrix*Satellite.solarpanels.D
            k = RMatrix*Satellite.solarpanels.k          
            t = SolarPanelAngle(j) 
            spanelN =  size(idxtop)(1,2)       
            qA =  CL_rot_defQuat(zeros(1,spanelN),A(:,idxtop)-D)
            qB =  CL_rot_defQuat(zeros(1,spanelN),B(:,idxtop)-D) 
            qC =  CL_rot_defQuat(zeros(1,spanelN),C(:,idxtop)-D) 
            qF =  CL_rot_defQuat(zeros(1,spanelN),F(:,idxtop)-D)
            q2 = CL_rot_defQuat(repmat(cos(t/2),1,spanelN),k*sin(t/2))      
            A(:,idxtop) = imag(q2*qA*q2')+D
            B(:,idxtop) = imag(q2*qB*q2')+D
            C(:,idxtop) = imag(q2*qC*q2')+D
            F(:,idxtop) = imag(q2*qF*q2')+D
            A(:,idxbot) = A(:,idxtop) 
            F(:,idxbot) = F(:,idxtop)
            B(:,idxbot) = C(:,idxtop) 
            C(:,idxbot) = B(:,idxtop)                                    
        end       
        x = []
        y = []
        z = []
        for i = 1 : Nplate
            x = [x,[A(1,i);B(1,i);F(1,i);C(1,i)]]
            y = [y,[A(2,i);B(2,i);F(2,i);C(2,i)]]
            z = [z,[A(3,i);B(3,i);F(3,i);C(3,i)]]
        end
        delete()
        plot3d(x,y,z)
        idGif = animaGIF(gcf(), idGif)
    end
    animaGIF(idGif)
endfunction

global Satellite
CL_init()



/*----------------------------------------------------------------------
Below aresome satellite configurations that I have coded in. To define a
Satellite, create a satellite object with the following variables

A,B,C: Coordinates of the A,B and C point of the plates of the satellite
Name: Satellite name
COM: Position of the center of mass
MOI: 3x3 moment of inertia matrix of the satellite
Mass: Mass of the satellite
surfacetemperature: Satellite surface temperature

Optional hardware that need to be set to %f if not using:
Magnetorquer: An object describing the magnetorquer hardware. The only
neccessary parameter in this object is the magnetorquer gain
solarpanels: An object describing the solarpanels hardware. The 
neccessary paramters are:
    D, E: points defining the rotation points and axis of the solarpanel
    align: option to align solarpanels, set to %t if you want alignment
    idxtop: index of top plate of the solarpanels
    idxbot: index of the bottom plate of the solarpanels    
----------------------------------------------------------------------*/

//shuttlecock configuration
A = [-5,0,0;0,1,0;0,0,1;-5,0,0;-5,1,0;-5,1,0;-5,0,1;-5,0,1;-5,0,0;-5,0,0;-5,1,0;-5,1,0;-10,6,1;-10,6,1]
C = [-5,0,1;0,0,0;-5,0,1;-5,1,0;-5,0,0;0,1,0;-10,0,6;-5,1,1;-10,-5,0;-5,0,1;-5,0,0;-10,1,-5;-5,1,1;-10,6,0]
B = [0,0,0;0,1,1;0,1,1;-5,0,1;0,1,0;-5,1,1;-5,1,1;-10,0,6;-5,0,1;-10,-5,0;-10,1,-5;-5,0,0;-10,6,0;-5,1,1]
Magnetorquer = struct("gain",1.e8)
Shuttlecock = struct("Name","Shuttlecock","COM",[-2.5,0.5,0.5],"MOI",10*eye(3,3),"Mass",100,"A",A,"B",B,"C",C,"surfacetemperature",280,"Magnetorquer",Magnetorquer,"solarpanels",%f)

//VELOX-C1
A = [0,0,0;0.75,0,0;0.75,0.50,0;0,0.50,0;0.75,0.50,-0.53;0.75,0.50,0;0,0,-0.53;0,0,-0.53;0,0.50,-0.53;0,0.50,-0.53]
B = [0,0,-0.53;0.75,0,-0.53;0.75,0.50,-0.53;0,0.50,-0.53;0.75,0,-0.53;0,0.50,0;0.75,0,-0.53;0,-0.50,-0.53;0,1,-0.53;0.75,0.50,-0.53]
C = [0.75,0,0;0.75,0.50,0;0,0.50,0;0,0,0;0,0.50,-0.53;0.75,0,0;0,-0.50,-0.53;0.75,0,-0.53;0.75,0.50,-0.53;0,1,-0.53]
D = [0.375,0,-0.53;0.375,1,-0.53]
E = [0.375,-0.50,-0.53;0.375,0.50,-0.53]
Magnetorquer = struct("gain",1.e6)
solarpanels = struct("D",D,"E",E,"align",%t,"idxtop",[7,9],"idxbot",[8,10])
MOI = [5.44,0,0;0,8.64,0;0,0,8.33]
Velox_C1 = struct("Name","Velox_C1","COM",[0.375,0.25,-0.265],"MOI",MOI,"Mass",123,"A",A,"B",B,"C",C,"Magnetorquer",Magnetorquer,"surfacetemperature",280,"solarpanels",solarpanels)

//VELOX-C1 modified
A = [0,0,0;0.75,0,0;0.75,0.50,0;0,0.50,0;0.75,0.50,-0.53;0.75,0.50,0;0,0,-0.265;0,0,-0.265;0,0.50,-0.265;0,0.50,-0.265]
B = [0,0,-0.53;0.75,0,-0.53;0.75,0.50,-0.53;0,0.50,-0.53;0.75,0,-0.53;0,0.50,0;0.75,0,-0.265;0,-0.50,-0.265;0,1,-0.265;0.75,0.50,-0.265]
C = [0.75,0,0;0.75,0.50,0;0,0.50,0;0,0,0;0,0.50,-0.53;0.75,0,0;0,-0.50,-0.265;0.75,0,-0.265;0.75,0.50,-0.265;0,1,-0.265]
D = [0.375,0,-0.265;0.375,1,-0.265]
E = [0.375,-0.50,-0.265;0.375,0.50,-0.265]
Magnetorquer = struct("gain",1.e6)
solarpanels = struct("D",D,"E",E,"align",%t,"idxtop",[7,9],"idxbot",[8,10])
MOI = [5.44,0,0;0,8.64,0;0,0,8.33]
Velox_C1mod = struct("Name","Velox_C1mod","COM",[0.375,0.25,-0.265],"MOI",MOI,"Mass",123,"A",A,"B",B,"C",C,"Magnetorquer",Magnetorquer,"surfacetemperature",280,"solarpanels",solarpanels)

//cubesat
A = [0,0,0;0,0,0;0,0.1,0;0.1,0,0;0,0,0.1;0,0,0]
B = [0.1,0,0;0,0.1,0;0.1,0.1,0;0.1,0,0.1;0,0.1,0.1;0,0,0.1]
C = [0,0.1,0;0,0,0.1;0,0.1,0.1;0.1,0.1,0;0.1,0,0.1;0.1,0,0]
Magnetorquer = struct("gain",1.e3)
Cubesat = struct("Name","Cubesat","COM",[0.05,0.05,0.05],"MOI",eye(3,3)/600,"Mass",1,"A",A,"B",B,"C",C,"Magnetorquer",%f,"surfacetemperature",280,"solarpanels",%f)

//flatplate
A = [0,0,0;0,0,0]
B = [1,0,0;0,1,0]
C = [0,1,0;1,0,0]
flatplate = struct("Name","flatplate","COM",[0.5,0.5,0],"MOI",eye(3,3)/600,"Mass",1,"A",A,"B",B,"C",C,"Magnetorquer",%f,"surfacetemperature",280,"solarpanels",%f)
