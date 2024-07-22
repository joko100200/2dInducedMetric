# 2dInducedMetric
 Calculates geodesics on 2d hypersurfaces. This code generates a uniform 2d flat mesh which is then implanted into 3d euclidean space based on three parmeterization functions. This code is written in fortran.
## How to run 
 Have gfortran installed
 type in terminal 
```bash
 ./manifold.sh
```
Which compiles the mesh creation part of the code. Which will output an executable runmesh.
```bash
 ./runmesh
```
 This generates the mesh that is based on your input parameters but it doesn't
 solve the geodesic equation.
 To run the geodesic equation you will need to compile
 ```bash
 ./Integrator.sh
 ```
 Which compiles the Integration part of the code.
 ```bash
 ./runIntegrator
 ```
 The above command runs the outputted executable.
 This will solve the geodesic equation given the parameters you put in the paramater.inp file.
 For a more indepth look into the variables specified by parameter.inp read the comments right
 below the variable declerations in either mesh.f95 or GeodesicIntegrator.f95.

## Visualising
 After executing runmesh you can you will have multiple output files.
 uv.out, xyzuv.out, connections.out, metric.out
 The quick explination of these files are: <br/> <br/>
   1. uv.out: Which is the 2d flat mesh with each row being a point where the two columns are (u,v) <br/>
   2. xyzuv.out: Which is the uv.out mesh once it has been through the parameterization defined in xyz.f95 his file is organized with each row being a single point and each colunm (x,y,z)  <br/>
 Now the two files below store components of the objects only at each point. The "tangent" (points) where these objects are valued at in the files is implicitely stored in the row of the data. More information about this [below](https://github.com/joko100200/2dInducedMetric?tab=readme-ov-file#how-the-mesh-works)  <br/>
   3. connections.out: This file stores the numerical values of connection which there are ofcourse in 2d are 6 unique components. Each row is a single point in space and the colunms are (111,112/121,122,211,212/221,222) <br/>
   4. metric.out This file sotres the numerical values of the metric and inverse metric at every point on the mesh where the colunms are (11,12/21,22,inv11,inv12/21,inv22)<br/><br/>

 
 
 # How the Mesh works
 
 Then to see the graph of both the 2d->3dmesh and the geodesic solution you will need to run gnuplot
 and the type in:
 set ticslevel 0
 splot 'xyzuv.out','xyzpath.out' with lines lw 3 
 This plots a 3d mesh and changes the line with to 3 pts
 Make sure that this is only for largeish values of hu and hv in the parameters file if not gnuplot will crash
 The way I suggest is to run the mesh first with hu hv = 0.1 or something large like that then run python3 addspacing.py
 when you do that it will create a lower res copy called 'visualxyzuv.out' and sometimes you will be able to use
 it with pm3d which will give it a color gradient based on the z value which is neat
 after that you will change the values of hu and hv to something much smaller like 0.005 run the mesh and then
 run the integrator. When you do this you will run gnuplot and type:
 set ticslevel 0
 splot 'visualxyzuv.out' with pm3d, 'xyzpath.out' with lines lw 3 lt rgb "green"
 This forces the xyzpath to be a neon green
 if you get an error or something looks messed up because of the mesh get rid of the pm3d part in the visual



Parameters file:
 Warnings rememeber this code forces all steps to be on the mesh. This means that if you
 put the velocity to small the point won't even move it must be able to make the gap.
 These are just warnings if you want to find out what each variable is please look at the code.
 There it will tell you what each one does unless they are a dummy variable in which case I might have been lazy and didn't put they were a dummy variable
 You can also refrence the refrences folder that includes a very indepth look into how the mesh works
 Even tho the mesh can be transformed into compact surfaces doesn't mean it works. Those require special boundry conditions
 which changes based on the transformation in xyz.f95 which I don't know how to do and I'm a bit to unintrested.
 So all surfaces you want to add must not wrap on itself i.e. They must not be compact

Code Explination:
 Better explination of the mesh can be found in refrences with the word file here I will just explain what code file does what

 mesh.f95: generates the mesh. ofc the mesh is an abstract and isn't real it just uses UVToIndex and IndexToUV to create all the nodes
 that mesh will evaluate the metric and connections on.
 Once mesh has finished it outputs the metric and the connections into there own files called metric.out and connections.out.
 In metric.out row i of the file is IndexToUV(i,u,v,...) is the uv point in the mesh and the 6 columns corespond to the 3 unique metric values and the 3 unique inverse metric
 1st | g_(11) g_(12) g_(22) g^(11) g^(12) g^(22)
 .
 .
 .
 Tst | g_(11) g_(12) g_(22) g^(11) g^(12) g^(22)
 The way the connections file is organized in the same way where the rows corespond to the u,v point by Index to U,V and we store all six unique connections as such.
 1st | 111 112 122 211 212 222 
 .
 .
 .
 Tth | 111 112 122 211 212 222
 
 T here is just used to represent the total number of points in the mesh which is calculated quite early on the mesh.
 If you open up mesh you will be able to see what I meah its ut*vt from this
 This code uses finite diffrence to calculate the partial derivatives
 
 xyz.f95: These are our functions where we take our flat 2d uv surface and we put all points through and write that output to a specified file. When the mesh code above calls it
 it gives it 'uv.out', T, 'xyzuv.out' xyz.f95 then takes uv.out reads every row u,v and then puts it through whatever function switch told it to then writes that output to xyzuv.out
 where xyzuv.out is organized by x(u,v),y(u,v),z(u,v) mesh.f95 reads that file after the call then puts all data into an array of size (T,3) to do finite diffrence derivatives.
 To see more explination go to the actual code it has explination of all the variables
 
 GeodesicIntegrator.f95: This code is what actually calculates the geodesic on the mesh because of the way mesh outputs the connections as well as the metric file we run this code seperatly which
 lets us quickly see results. Basically we don't need to recreate the mesh everytime we want to run the geodesic calculation we can just run the geodesic calculation.
 This code is completely intrinsic and the only refrence to xyz.f95 is to make uvpath.out -> xyzuvpath.out which is how we can see the curvature. However you can graph uv.path by itself and see the curve
 in 2d. The way it integrates the geodesic equation is as follows:
 Reads initial conditions for initial u and v
 Reads initial conditions for initial u and v velocities
 then plugs it into the Geodesic equation to solve for u and v acceleration (you can find this in GeodesicIntegrator it will be in the main integration loop with all the connections(i,1) etc)
 then velocity u (new) = velocity u (old) + acceleration u (just calculated) * dt (from input)
 same for v
 then position u (new) = position u (old) + velocity u (new) * dt
 same for v
 then we force these new v and u points to be on the mesh by putting it through our two mesh subroutines
 UVToIndex()
 IndexToUV()
 then we write this to uv.path file write("uvpath.out",*) u,v
 then u (old) = u (new)
 same for v
 and the loop contines for time (input parameter) steps
 
 Once the main loop is over with we create two slight diffrences in the path using a function that is zero at the first step and zero at the last step.
 You can see these functions and the main one once you finish running all the codes and opening gnuplot and typing:
 plot 'uvpath.out' with lines, 'uv+path.out' with lines, 'uv-path.out' with lines
 We then put these 3 uv paths through a distance subroutine (I did it with a function and it was jank so I will be using subroutines from now on). Ofc the path with the smallest distance must be
 our geodesic path given in 'uvpath.out'. It will write to screen when Integrator is finished distance('uv+path.out') distance('uvpath.out') distance('uv-path.out'). The way the distance calculator works
 Is just a simple line integral you can see it all the way at the bottom of GeodesicIntegrator inside the subourtine distance S = S + sqrt(...) all inside one bigger loop which is the integration loop.

