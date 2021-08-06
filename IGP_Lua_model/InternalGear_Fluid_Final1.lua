---------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------Internal Gear Pump--------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------

-- This script is intended to generate the Parametrical components of the internal gear pump.
-- Three main components of the Internal Gear Pump  the Internal Gear, External Gear / Pinion,
-- along with the cresant filler and their dependiencies are modelled. This generated model 
-- provides the user  with the flexibility to explore the wide change in parameters and their 
-- corresponding interdependiencies. Further seperate files for each compoments are provided 
-- along with the generated '.stl' files and 'Gcodes' which could be help incorporating models
-- in further work or could be used for 3d printing of the same.  
---------------------------------------------------------------------------------------------------------------------------------
--COMMENT LINE 443 (Front Flange) TO VIEW THE COMPONENTS OF INTERNAL GEAR PUMP




---------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------Declaration and Calculation End--------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------

enable_variable_cache = true;                                                                                               -- caches the pointer

----------------------------------------------Passing parameters to user interface-----------------------------------------------
z_n=ui_numberBox("Number Of Teeth",25); 						                                                            -- Input for number of teeth, ideal: from 25
m=ui_numberBox("Module Of Gear",3);								                                                            -- Module of gear, Ideal Input: 3
alpha_t=ui_scalarBox("Pressure Angle(Deg)",20,1);				                                                            -- Pressure angle (changes the meshing parameters)
width=ui_numberBox("Width(mm)",15);								                                                            -- Width/Thickness of the gear
x_coef_int=ui_scalarBox("Internal Profile Shift(mm)",-0.2,0.1);	                                                            -- Profile shift factor for internal gear
x_coef_ext=ui_scalarBox("External Profile Shift(mm)",0,0.1);	                                                            -- Profile shift factor for external gear
h_a_coef_p=ui_scalarBox("Addendum Coefficient(mm)",1,0);					                                                -- Addendum height factor
h_f_coef_p=ui_scalarBox("Dedendum Coefficient(mm)",1.25,0);					                                                -- Dedendum height factor
rotation = ui_numberBox("Fluid Flow Rotation",0)*0.1;                                                               			-- For rotation
transform_pick = ui_pick('Drag the top flange to move it')


--------------------------Function Ellipse--------------------------------------------------------------------------------------
function ellipse(a,b,n)
   -- example for making a contour that will be used with linear_extrude to build a shape
   -- the contour is ellipse with the half-axis a,b build out of n support points
   -- the contour is a table of vectors (IceSl) having n+1 elements start points = end point!
   local XY= {}                                   -- create the matrix => table, use local for local variables!
   for i= 1,n do                                  -- loop matlab for i=1:n
      XY[i]= v(math.cos(2*math.pi*(i-1)/n)*a,     -- using IceSl vector command "v" to fill table with vectors
               math.sin(2*math.pi*(i-1)/n)*b,5)   -- the contour is not closed: start point and end point differs
   end                                            -- end loop  
   XY[n+1]= XY[1]                                 -- close contour set end point to start point
   return XY                                      -- return value of function: table of vectors
end  


--------------------------Function definition for Parametric equations for the involute profile points---------------------------
 

function tooth_involute(base_radius, inv_alpha)                                                                             -- inv_alpha = Involute angle
    return v(base_radius*(math.sin(inv_alpha) - inv_alpha*math.cos(inv_alpha)), 
			base_radius*(math.cos(inv_alpha) + inv_alpha* math.sin(inv_alpha)))
end

-----------------------------------Function defining for mirroring or inverting profile points-----------------------------------
-- coord -Co-ordinates
function tooth_mirror(coord) 
	return v(-coord.x, coord.y) 
end

---------------------------------------------Function defining the rotational matrix---------------------------------------------

function rotate_points(angle, coord)                                                                                        -- angle = angle, coord = co-ordinates
    return v(math.cos(angle) * coord.x + math.sin(angle) * coord.y, math.cos(angle) * coord.y - math.sin(angle) * coord.x)
end

--------------------------------------Function defining angle between corresponding points---------------------------------------

function involute_angle(r_p1, r_p2)                                                                                       
	return (math.sqrt((r_p2 * r_p2 - r_p1 * r_p1) / (r_p1 * r_p1))) 
end

----------------------------Function Working Pressure angle -------------------
function wkp(x)
	return (((math.pow(3,(1/3)))*(math.pow(x,(1/3)))) - (2*x/5) + (math.pow(9/175*3,(2/3)))*(math.pow(x,(5/3))) - (math.pow(2/175*3, (1/3)))
						*(math.pow(x,(7/3))) - ((144/67375)*(math.pow(x,(9/3))) + (3258/3128215)*(math.pow(3,(2/3)))*(math.pow(x,(11/3)))
						- (49711/153278125)*(math.pow(3,(1/3)))*(math.pow(x,(13/3))) - (1130112/9306171875)*(math.pow(x,(15/3)))
						+ (5169659643/95304506171875)*(math.pow(3,(2/3)))*(math.pow(x,(17/3)))))
end



------------------------Function for the calculation for the required Internal and External Gear Profiles------------------------

function gear_profile(z, m, alpha_t, x_coef, h_a_coef, h_f_coef, width)

    local inv_xy = {}
	
-- Definition of the input parameters and calculation of the other base parameters for the involute profile.
	
    m_t = m;  						                                                                                        -- Module of gear

    alpha_t_rad = alpha_t*math.pi/180;                                                                                      -- Pressure angle
     
    z_t = z;  						                                                                                        -- Number of teeth

    x_coef = x_coef;				                                                                                        -- Profile shift co-efficient

    --x_coef_ext = 0; 				                                                                                        -- Profile shif coeff external gear

    c = 0.167 * m_t;  			                                                                                            -- Clearance of tooth
	--c = 0.5 * m_t; 
	h_acoef = h_a_coef; 	                                                                                                -- Addendum coefficient
	
	h_dcoef = h_f_coef; 		                                                                                            -- Dedendum coefficient

-- Base formulae calculations required for the gear profile

-- Pitch Diameter:
    d_p = z_t * m_t 	   			                                                                                        -- Pitch Diamter
    r_p = d_p / 2 				   	                                                                                        -- Pitch radius

    
-- Base Diameter:
    d_b = d_p * math.cos(alpha_t_rad);                                                                        -- Base diameter of gear
	--d_b = d_p * math.cos(alpha_t_rad) + 2*x_coef*m_t                                                                        -- Base diameter of gear
    r_b = d_b / 2 				                                                                                            -- Base radius 


-- Addendum and Dedendum
	h_a = m_t* h_acoef; 			                                                                                        -- Addendum

	h_f = m_t* h_dcoef  			                                                                                        -- Dedendum

-- Tip Diameter / Addendum Diameter 
	
	d_a = m_t*(z_t+2*x_coef+2);
	--d_a = (d_p + (2 * m_t * (1 + x_coef)))                                                                                  -- Addendum diameter

    r_a = d_a / 2 				                                                                                            -- Addendum radius
  
-- Root Diameter / Dedendum Diameter
	d_f = m_t*z_t+ 2*x_coef*m_t - 2*(h_f)                                                                                   -- Root diamter
    r_f = d_f / 2 			                                                                                                -- Root_radius

-- True Involute Diameter:

    Q_fp = ((m_t*((math.pi/4)-(math.tan(alpha_t_rad)))-c*math.tan(alpha_t_rad))*(1+math.sin(alpha_t_rad)))/(math.cos(alpha_t_rad))
	d_TIF = math.sqrt(math.pow(d_p * math.sin(alpha_t_rad) - 2 *(h_a - (m_t * x_coef) - Q_fp *(1 - math.sin(alpha_t_rad)))/(math.sin(alpha_t_rad)), 2) +   d_b * d_b)
	--d_TIF = d_TIF-2*c
	--d_TIF = math.sqrt(math.pow(d_p * math.sin(alpha_t_rad) - 2 *(h_a - (m_t * x_coef) - h_f *(1 - math.sin(alpha_t_rad))), 2) +   d_b * d_b);                                -- True involute diameter
    r_TIF = d_TIF / 2; 				    -- true form radius 

--- Tooth thickness calculation on the form circle:
	
-- Tooth thickness on the pitch Circle

    S_0 = m_t*((math.pi/2)+ 2*x_coef* math.tan(alpha_t_rad));

-- Involute function 

    inv_a = math.tan(alpha_t_rad) - alpha_t_rad; 
	
-- Tooth thickness on the form circle 

    alpha_f = math.acos((d_p*math.cos(alpha_t_rad))/d_TIF); 				                                                -- Involute angle along form circle
    inv_Ff  = math.tan(alpha_f) - alpha_f;									                                                -- Involute function
    s_f     = d_TIF*((S_0/d_p)+inv_a - inv_Ff)+(x_coef*math.tan(alpha_f));                                                  -- Thickness of the gear teeth along form circle 
    omega   = (s_f)/(0.5*d_TIF); 											                                                -- Angle swept along form circle for corresponding thickeness
	
-- Function for stariing and ending of involute between two radius

    tooth_ang = (((math.pi * m_t / 2) + 2 * m_t * x_coef * math.tan(alpha_t_rad)) /
                    r_p + 2 * math.tan(alpha_t_rad) - 2 * alpha_t_rad)
    

    res = 30;                                                                                                               -- Defining iterating points 
 
-- Iterating the points for the involute points and iterating teeth around the circle

    for i = 1, z_t do
        th = 2 * math.pi 				                                                                                    -- Iteration angle of the teeth around the gear diameter
        th1 = involute_angle(r_b, r_TIF);                                                                                     -- Start angle to define the start of the involute curve 
        th2 = involute_angle(r_TIF, r_a);                                                                                     -- End asngle to determine the end of the involute curve 

-- Defining iterations for the involute teeth profile 

        for j = 1, res do

            --inv_xy[#inv_xy + 1] = rotate_points(th*i/z_t, tooth_involute(r_b, (th1 +(th2-th1) * j / res)))				    -- The one side of the involute points are obtained 
			inv_xy[#inv_xy + 1] = rotate_points(th*i/z_t, tooth_involute(r_TIF, ((th2) * j / res)))
        end
        for j = res, 1, -1 do

            inv_xy[#inv_xy + 1] = 
			--rotate_points(th*i/z_t, rotate_points(omega, tooth_mirror(tooth_involute(r_b, (th1 + (th2 - th1)* j / res)))))  -- The second side of the involute points are obtained using mirroring function 
			rotate_points(th*i/z_t, rotate_points(omega, tooth_mirror(tooth_involute(r_TIF, ((th2)* j / res)))))
        end

    end

    inv_xy[#inv_xy + 1] = inv_xy[1] 	                                                                                    -- Used to generate a table of values 
	interior = linear_extrude(v(0,0,width),inv_xy)
    return interior
end
------------------------------Defining the Default Parameters for the Internal Gear Pump Components------------------------------
function gear(params)
      return gear_profile(
			params.z or 20,				                                                                                    -- Number of teeth of the Internal Gear
			params.m or 3, 			                                                                                        -- Module  
			params.alpha_t or 28, 		                                                                                    -- Pressure Angle
            params.x_coef or 0, 		                                                                                    -- Profile shift coefficient 
			params.h_a_coef or 1,		                                                                                    -- Addendum coefficient
			params.h_f_coef or 1.25,                                                                                        -- Dedendum coefficient
			params.width or 15)			                                                                                    -- Thickness/Width 
end

-----------------------------------Function defining the parametrical equation for the circle------------------------------------
-- r -radius of circle
function circle(r)
  local x, y = 0, 0
  local XY={}
  for i = 1, 360 do
    local angle = i * math.pi / 180
    XY[i] = v(x + r * math.cos( angle ), y + r * math.sin( angle ))
  
  end  
  return XY
end 

----------------------------------------------------For Internal Gear Extrude----------------------------------------------------

function extrude(Contour, angle, dir_v, scale_v, z_steps)
-- extrude a Contour to a shape by turning the contour to angle in z_steps
-- extrude a Contour in a dircetion given by the vector dir_v
-- extrude a Contour with a scaling factor given by vector scale_v 
-- Contour: a table of vectors as a closed contour (start point and end point is the same)
-- angle: roation angle of contour along z_steps in deg
-- dir_v: vector(x,y,z) direction of extrusion
-- sacle_v: vector(x,y,z) scaling factors for extrudion
-- z_steps: number of steps for the shape, mostly z_steps=2 if angle is equal zero
 
   local n_cont= #Contour
   local angle= angle/180*math.pi
   local Vertex= {}
   
   for j= 0,z_steps-1 do
      local phi= angle*j/(z_steps-1)
      local dir_vh= dir_v*j/(z_steps-1)
      local scale_vh= (scale_v - v(1,1,1))*(j/(z_steps-1)) + v(1,1,1)
	  for i= 1,n_cont-1 do
          Vertex[i+j*n_cont]= v((dir_vh.x + scale_vh.x * (Contour[i].x*math.cos(phi) - Contour[i].y*math.sin(phi))),
                                (dir_vh.y + scale_vh.y * (Contour[i].x*math.sin(phi) + Contour[i].y*math.cos(phi))),
                                (dir_vh.z * scale_vh.z))
      end
	  
      table.insert(Vertex,Vertex[1+j*n_cont])
   end

   local vertex_sum_1 = v(0,0,0)
   local vertex_sum_m = v(0,0,0)
   
   for i= 1,n_cont-1 do
      vertex_sum_1= vertex_sum_1 + Vertex[i]
      vertex_sum_m= vertex_sum_m + Vertex[i+n_cont*(z_steps-1)]
   end    

   table.insert(Vertex,vertex_sum_1/(n_cont-1))                                                                             --n_cont*m_cont + 1
   table.insert(Vertex,vertex_sum_m/(n_cont-1))                                                                             --n_cont*m_cont + 2


   Tri= {}                                                                                                                  -- !!! the index on table with Vertex starts with zero !!!
   local k= 1
   
   for j=0,z_steps-2 do
      for i= 0,n_cont-2 do
         Tri[k]=   v(i, i+1, i+n_cont) + v(1,1,1)*n_cont*j
         Tri[k+1]= v(i+1, i+n_cont+1, i+n_cont) + v(1,1,1)*n_cont*j
         k= k+2
      end
   end
   for i= 0,n_cont-2 do
      Tri[k]= v(i+1,i,n_cont*z_steps)
      k= k+1
   end
   for i= 0,n_cont-2 do
      Tri[k]= v(i+n_cont*(z_steps-1),i+1+n_cont*(z_steps-1),n_cont*z_steps+1)
      k= k+1
   end
   return(polyhedron(Vertex,Tri))
end

test = x_coef_int - x_coef_ext;
if test < -0.1 then
	x_coef_ext = x_coef_int + 0.1;
-- elseif x_coef_int < -0.3 then
	-- x_c = x_coef_int + 0.1;
	-- x_coef_ext = x_c - 0.1;
-- else
	-- x_coef_ext = x_coef_ext;
end	
if alpha_t < 4 then
	alpha_t = 5;
end 
----------------------------------------------Function for Internal Gear Formation-----------------------------------------------
function internal_gear()

  internalGear = gear({z=z_n;m=m;alpha_t=alpha_t;x_coef=x_coef_int;h_a_coef=h_a_coef_p;h_f_coef=h_f_coef_p;width=width})

end
internal_gear()



---------------------Calculation for the centre distance between gears using the profile shift coefficients----------------------  
   
-- Involute function -Operating pressure angle
x_coef_int = x_coef_int; 		                                                                                            -- Profile shift coefficient internal gear 
x_coef_ext = x_coef_ext; 		                                                                                            -- Profile shift coefficient External gear 
z_2 = z_n;                 		                                                                                            -- Number of teeth Internal Gear
z_1 = z_n-8;				 		                                                                                        -- Number of teeth External Gear

alpha_rad = alpha_t*math.pi/180                                                                                             -- Preassure angle 
inv_a = math.tan(alpha_rad) - alpha_rad;                                                                                    -- Involute function 
inv_aw = ((2*math.tan(alpha_rad) * (x_coef_int - x_coef_ext))/(z_2 - z_1)) + inv_a;                                         -- Involute function working pressure angle 

-- Working preassure angle
--alpha_aw = alpha_awn*math.pi/180;                                                                                           -- Working pressure angle 
alpha_aw = wkp(inv_aw)
-- Centre distance incremental factor 
y_c = ((z_2 - z_1) * (math.cos(alpha_rad) - math.cos(alpha_aw)))/ (2* math.cos(alpha_aw));

-- Center distance between Gears
a_x = (((z_2 - z_1)/2)+y_c) * m; 	                                                                                        -- Center distance 

meshDistance = a_x;				                                                                                            -- Center distance is assagined for cresant calculation
---------------------------------------------------------------------------------------------------------------------------------

addOn_distance= m * 2;
crescent_outerRadius=r_TIF;		                                                                                            -- Radius of the internal gear for the cylinder 
crescent_clearanceOut=c;		                                                                                            -- Clearance of the internal gear for the cylinder
crescent_rootRadius=r_f;		                                                                                            -- Root radius of the internal gear
internal_radius=r_a;																										-- Pitch radius of internal gear
crescent_cylTop=ccylinder(crescent_outerRadius-crescent_clearanceOut*2,width);	                                            -- Top cylinder for crescent 
radius=(z_n*m)/2;

----------------------------------------------Function for External Gear Formation-----------------------------------------------
function external_gear()
  
  externalGear = gear({z=z_n-8;m=m;alpha_t=alpha_t;x_coef=x_coef_ext;h_a_coef=h_a_coef_p;h_f_coef=h_f_coef_p;width=width})
  
end
external_gear()

--------------------------------------------------------Crescent Filler----------------------------------------------------------

crescent_innerRadius=r_a;   			                                                                                    -- Radius of the external gear for the cylinder
crescent_clearanceIn=c;					                                                                                    -- Clearance of the external gear for the cylinder
 
crescent_translation=meshDistance;		                                                                                    -- Translation for crescent formation

crescent_cylBottom=translate(0,crescent_translation,0)*ccylinder(crescent_innerRadius+crescent_clearanceIn*2,width); 	    --Inner circle for internal gear formation

crescent_cubeRight=translate(crescent_innerRadius+crescent_clearanceOut,meshDistance+m*2,0)*
                   cube(crescent_rootRadius+crescent_clearanceOut, crescent_rootRadius+crescent_clearanceOut,width+0.1)     -- Cube to remove sharp edges of the crescent
crescent_cubeLeft=translate(-(crescent_innerRadius+crescent_clearanceOut),meshDistance+m*2,0)*
                   cube(crescent_rootRadius+crescent_clearanceOut, crescent_rootRadius+crescent_clearanceOut,width+0.1)     -- Cube to remove sharp edges of the crescent
crescent_Main=translate(0,0,width/2+0.1)*difference(crescent_cylTop,crescent_cylBottom)                                     -- Crescent Filler Formation with sharp edges
crescent_Full=intersection(difference(crescent_Main,crescent_cubeRight),
                   difference(crescent_Main,crescent_cubeLeft))                                                             -- Crescent Filler Formation without sharp edges 
				   
				   
 
---------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------Declaration and Calculation End--------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------
r1 = rotate(0,0, rotation)                                                                                                  -- Internal Gear Rotation
r2 = rotate(0,0, rotation*z_2/z_1)                                                                                          -- External Gear Rotation

---------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------Emitting Section---------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------

---------------------------------------------------Emitting the Internal Gear----------------------------------------------------
radius=((z_n)*m/2);                                                                                                         -- Radius for outer cylinder
outer_circle= extrude(circle(internal_radius-0.1), 0, v(0,0,width), v(1,1,1), 20)                                       		-- Outer circle for internal gear formation

emit(r1*difference(outer_circle,internalGear),128);                                                                         -- Internal gear formation
emit(translate(0,0,-0.1)*cylinder(internal_radius,2),128)																		-- Internal gear Base formation
intr=r1*difference(outer_circle,internalGear)
---------------------------------------------------Emitting the External Gear and shaft----------------------------------------------------

bore_externalGear= extrude(circle(addOn_distance*2), 0, v(0,0,width), v(1,1,1), 20)	        								-- Bore formation of external gear
flangeCube=translate(0,-meshDistance,0)*cube(addOn_distance*2,addOn_distance*2,width)
flangeCylinder=cylinder(addOn_distance*2,width+50)
cube_flange=translate(0,-meshDistance,0)*cube(addOn_distance*2,addOn_distance*2,width+5)
flangeDiff=difference(bore_externalGear,flangeCube)
flangeFull=difference(externalGear,flangeDiff)
fln=translate(0,meshDistance,0.1)*r2*difference(flangeCylinder,cube_flange)
emit(translate(0,meshDistance,0.1)*r2*difference(flangeCylinder,cube_flange),80)										 	-- Shaft Formation with key-slot
ex=translate(0,meshDistance,0)*r2*difference(flangeFull,flangeDiff)
emit(translate(0,meshDistance,0)*r2*difference(flangeFull,flangeDiff),2)                         						 	-- External Gear Formation

--------------------------------------------------Emitting the Crescent Filler---------------------------------------------------
cr=crescent_Full
emit(crescent_Full)                                                                                                        	-- Crescent Filler Formation

-------------------------------------------------------Casing------------------------------------------------------

back_flange=cylinder(internal_radius+16,0.5)				                                                            
casing_cyl=difference(cylinder(internal_radius+16,width),cylinder(internal_radius+1.75,width))												-- Casing for gear					
																					
case_upper=rotate(45,Z)*translate(0,0,width+0.2)*cylinder(internal_radius+16,0.1)																-- Bore for ports in the casing
hole_casing=translate(0,meshDistance,width+0.2)*extrude(circle(addOn_distance*2+0.5), 0, v(0,0,0.1), v(1,1,1), 20)	        -- Bore for ports in the casing

-----------------------------------------------------------Fluid--------------------------------------------------------------------
ayy=extrude(circle(internal_radius+6), 0, v(0,0,width), v(1,1,1), 20)	
--a=(cube(2*internal_radius+10,2*internal_radius+10,width)) 
--emit(a)   	
--r1*difference(difference(r2*difference(difference(a,crescent_Full),ex),intr),fln	
gg=difference(difference(ayy,intr),ex)
kuy=difference(difference(gg,fln),cr)
--emit(gg,5)													
--emit(difference(difference(gg,fln),cr),50)
intdiff=difference(cylinder(internal_radius,width),cylinder(crescent_rootRadius,width))
yuk=difference(intdiff,intr)
--emit(yuk,80)
uik=difference(translate(0,meshDistance,0)*cylinder(crescent_innerRadius,width),ex)
--emit(uik,50)
----------------------------------------------diff colour fluids----------------------------------------------------------------------

trns=translate(rotation,0,0)
	cyl=rotate(270,X)*difference(cylinder(width/2+2,internal_radius*2), translate(0,-(width/2)+2,0)*cylinder(width/2-1,internal_radius*2))	--piping
	emit(rotate(45,Z)*translate(0,internal_radius+2+31.5,3)*cyl,22)
	cyl1=rotate(270,-Y)*difference(cylinder(width/2+2,internal_radius*2), translate(-(width/2)+2,0,0)*cylinder(width/2-1,internal_radius*2))	--piping
emit(rotate(45,Z)*translate(internal_radius+2+31.5,0,3)*cyl1,22)
if(rotation~=0) then


	

	
	
				
		-- end
	if(rotation>0) then
		thk5=translate(-internal_radius*2+2,0,0)*trns*cube(internal_radius*2+10,2*internal_radius+32, width)
		aj5 =(difference(thk5,cylinder(internal_radius+2,width-2)))
		ykj5= (difference(thk5,aj5))
		uk5=intersection(ykj5,kuy)
		emit(uk5,10)
	
		rk=internal_radius*2
		--inner_fluid1= rotate(270,X)*cylinder(width/2-2,30+internal_radius*2-rotation)
		inner_fluid1= rotate(270,X)*extrude(ellipse(width/2-4,width/2-2,100), 720, v(0,0,30+internal_radius*2-rotation), v(1,1,1), 200)
		emit(rotate(45,Z)*translate(0,internal_radius+2,width/2-1)*inner_fluid1,10)												-- left-side fluid inner flow fluid
		
		inner_fluid2= rotate(270,X)*extrude(ellipse(width/2-4,width/2-2,100), 1000, v(0,0,30+internal_radius*2-rotation), v(1,1,1), 200)
		emit(rotate(46,Z)*translate(1,internal_radius+2,width/2-1)*inner_fluid2,10)
		
		if(rotation>internal_radius*2 and rotation<=(internal_radius+2)*4) then
		
	
		inner_fluid3= rotate(270,-Y)*extrude(ellipse(width/2-4,width/2-2,100), 720, v(0,0,30+rotation-internal_radius*2), v(1,1,1), 200)											
		emit(rotate(45,Z)*translate(internal_radius+2,0,width/2-1)*inner_fluid3,10)												--right-side fluid outer flow
		
		inner_fluid4= rotate(270,-Y)*extrude(ellipse(width/2-4,width/2-2,100), 1000, v(0,0,30+rotation-internal_radius*2), v(1,1,1), 200)											
		emit(rotate(46,Z)*translate(internal_radius+2,-0.5,width/2-1)*inner_fluid4,10)
		
		end

	

	elseif(rotation<0) then
		
		thk6=translate(internal_radius*2+2,0,0)*trns*cube(internal_radius*2+10,2*internal_radius+32, width)
		aj6 =(difference(thk6,cylinder(internal_radius+2,width-2)))
		ykj6= (difference(thk6,aj6))
		uk6=intersection(ykj6,kuy)
		emit(uk6,10)
		
		
		inner_fluid_c= rotate(270,-Y)*extrude(ellipse(width/2-4,width/2-2,100), 720, v(0,0,30+internal_radius*2+rotation), v(1,1,1), 200)	
		emit(rotate(45,Z)*translate(internal_radius+2,0,width/2-1)*inner_fluid_c,10)												-- right-side fluid inner flow
		
		inner_fluid_c1= rotate(270,-Y)*extrude(ellipse(width/2-4,width/2-2,100), 1000, v(0,0,30+internal_radius*2+rotation), v(1,1,1), 200)	
		emit(rotate(45,Z)*translate(internal_radius+2,-0.5,width/2-1)*inner_fluid_c1,10)
	
		if(rotation<-(internal_radius-2)*2 and rotation>=-(internal_radius+2)*4) then
		
		
		
		outer_fluid1= rotate(270,X)*extrude(ellipse(width/2-4,width/2-2,100), 720, v(0,0,30-rotation-internal_radius*2), v(1,1,1), 200)	
		--outer_fluid1= rotate(270,X)*cylinder(width/2-2,30-rotation-internal_radius*2)
		emit(rotate(45,Z)*translate(0,internal_radius+2,width/2-1)*outer_fluid1,10)												--left-side fluid outer flow
		
		
		outer_fluid2= rotate(270,X)*extrude(ellipse(width/2-4,width/2-2,100), 1000, v(0,0,30-rotation-internal_radius*2), v(1,1,1), 200)	
		--outer_fluid1= rotate(270,X)*cylinder(width/2-2,30-rotation-internal_radius*2)
		emit(rotate(45,Z)*translate(0,internal_radius+2,width/2-1)*outer_fluid2,10)	
		end

	end

end




---------------------------------------------------------- Inlet port ----------------------------------------------------------
bore_inner=rotate(45,Z)*translate(radius+addOn_distance+16,0,width/2+0.2)*(rotate(270,Y)*cylinder(width/2-1,radius+addOn_distance+16))
bore_diffInner=difference(casing_cyl,bore_inner)
inner_port=difference(rotate(270,-Y)*cylinder(width/2-1,25),rotate(270,-Y)*cylinder(width/2-1.5,25))
inner_coupler= difference(rotate(270,-Y)*cylinder(width+2,3),rotate(270,-Y)*cylinder(width/2-1,3))
emit(rotate(45,Z)*translate(internal_radius+2,0,width/2)*inner_port,50)																		
emit(rotate(45,Z)*translate(internal_radius+2+25,0,width/2)*inner_coupler,50)
emit(rotate(45,Z)*translate(internal_radius+2+28.5,0,width/2)*inner_coupler,50)			
--------------------------------------------------------- Outlet port -----------------------------------------------------------
bore_outer=rotate(45,Z)*translate(0,radius+addOn_distance+16,width/2+0.2)*(rotate(270,-X)*cylinder(width/2-1,radius+addOn_distance+16))
bore_diffOuter=difference(casing_cyl,bore_outer)
output_bore=intersection(bore_diffOuter,bore_diffInner)
outer_port=difference(rotate(270,X)*cylinder(width/2-1,25),rotate(270,X)*cylinder(width/2-1.5,25))
outer_coupler= difference(rotate(270,X)*cylinder(width+2,3),rotate(270,X)*cylinder(width/2-1,3))
emit(rotate(45,Z)*translate(0,internal_radius+2,width/2)*outer_port,50)
emit(rotate(45,Z)*translate(0,internal_radius+2+25,width/2)*outer_coupler,50)																
emit(rotate(45,Z)*translate(0,internal_radius+2+28.5,width/2)*outer_coupler,50)
----------------------------------------------------------For screws--------------------------------------------------------------
hbase = 0.1
baseCylinder =translate(0,0,0)*case_upper

radiusHole = 2.5
heightHole = width+5
notcylinder = cylinder(radiusHole,heightHole)
totalHole = (translate(internal_radius+9,0,0)*notcylinder)

for i = 360/6 , 360 , 360/6 do
		totalHole = union{totalHole,rotate(0,0,i)*totalHole}
end


output_upperCase=(difference(baseCylinder,totalHole))
output_upperCase1=(difference(case_upper,hole_casing))															-- Front flange with bore for shaft
--emit(transform_pick*intersection(output_upperCase1,output_upperCase),50)														-- Front flange with screw holes									


radiusHole1 = 2.5
heightHole1 = width+5
holeCyl1 = cylinder(radiusHole1,heightHole1)
holeTotal1 = (translate(internal_radius+9,0,0)*holeCyl1)

for i = 360/6 , 360 , 360/6 do
		holeTotal1 = union{holeTotal1,rotate(0,0,i)*holeTotal1}
end
emit(difference(output_bore,holeTotal1),50)																		-- Gear casing with screw holes
			
radiusHole2 = 2.5
heightHole2 = width-5
holeCyl2 = cylinder(radiusHole2,heightHole2)
holeTotal2= (translate(internal_radius+9,0,0)*holeCyl2)

for i = 360/6 , 360 , 360/6 do
		holeTotal2 = union{holeTotal2,rotate(0,0,i)*holeTotal2}
end
emit(translate(0,0,-0.2)*difference(back_flange,holeTotal2),50)																		-- Back flange with screw holes	
--emit(translate(0,0,-10)*cube(3*radius+addOn_distance,4*radius+addOn_distance,0.1),60)                              -- Rectangular Base Plate Formation
--------------------------------------------------Setting Brushing for Colours---------------------------------------------------

set_brush_color(128,0,0,0)                                                                                                  -- Internal Gear
set_brush_color(100,0.7,0.7,0.1)                                                                                            -- External Gear
set_brush_color(80,0.29,0.26,0.24)                                                                                 			-- Flange
--set_brush_color(70,0.3,0.3,0.88)                                                                                            -- Base Circular Plate
set_brush_color(50,0.05,0.18,0.38) 
set_brush_color(55,0.308,0.166,0.104)																							-- Port
set_brush_color(60,0.1,0.1,0.11)                                                                                            -- Rectangular Base Plate
set_brush_color(30,0.3,0.6,0.9)

set_brush_color(10,0.300,0.600,0.900)
set_brush_color(15,0.60,0.947,0.942)
set_brush_color(8,0.315,0.535,0.963)
---------------------------------------------------Fluid Colour-----------------------------------------------------------------
--set_brush_color(21,0.167,0.172,0.180)
set_brush_color(21,0.455,0.599,0.838)
set_brush_color(22,0.422,0.449,0.493)
set_brush_color(23,0.192,0.365,0.631)
set_brush_color(20,0.89,0.17,0.89)

--set_brush_color(2,0.934,0.635,0.113)
set_brush_color(2,0.955,0.623,0.044)
----------------------------------------------------------Font Printing----------------------------------------------------------

-- f_n = font(Path .. 'font.ttf')
-- output = translate(0,-2*radius+addOn_distance-0.5,-3.1)*scale(0.3,0.3,1)*f_n:str('Optimal Gear Solutions',3)                     -- Font Location and Definition
-- emit(output)

--screenshot()
---------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------Emitting Section End-------------------------------------------------------