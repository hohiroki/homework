module variable

implicit none

integer:: coordinate_mode                    !1-Cartesian coordinate,2-Cylindrical cooridnate
real:: y_R0                                  !initial radius in the  y direction 
real:: x_length,y_length                     !length in the x and y direction
real:: rou,Cp                                !mass denstity and heat capacity at constant pressure

real:: k_boundary                            !thermal conductivity
real:: T_f                                   !Temperature of free fluid
real:: h                                     !convection coefficient
real:: q_boundary                            !heat flux boundary

integer:: i,j,step
integer,parameter:: x_node_number=101       !node number in x direction
integer,parameter:: y_node_number=101       !node number in y direction

real,dimension(2:x_node_number,2:y_node_number):: k !thermal conductivity

real,dimension(1:x_node_number):: x_node            !node dimension
real,dimension(2:x_node_number):: x_face            !control face position
real,dimension(2:x_node_number):: dx_node           !distance between the two adjacent nodes
real,dimension(2:(x_node_number-1)):: dx_face       !distance between the two adjacent control faces
real,dimension(2:x_node_number):: dx_face_node_E    !distance between the control face and its east node
real,dimension(2:x_node_number):: dx_face_node_W    !distance between the control face and its west node

real,dimension(1:y_node_number):: y_node            !node position
real,dimension(2:y_node_number):: y_face            !control face position
real,dimension(2:y_node_number):: dy_node           !distance between the two adjacent nodes
real,dimension(2:(y_node_number-1)):: dy_face       !distance between the two adjacent control faces
real,dimension(2:y_node_number):: dy_face_node_S    !distance between the control face and its south node
real,dimension(2:y_node_number):: dy_face_node_N    !distance between the control face and its north node
real,dimension(1:y_node_number):: y_R_node      !radius at the node in the y direction
real,dimension(2:y_node_number):: y_R_face      !radius at the control face in the y direction

real,dimension(1:x_node_number,1:y_node_number):: T1,T2                              !T1:former,T2:current
real,dimension(2:(x_node_number-1),2:(y_node_number-1)):: a_E,a_W,a_S,a_N,a_P,a_P0,b !influcence coefficient
real,dimension(2:(x_node_number-1),2:(y_node_number-1)):: S_c,S_p                    !coefficient for source term
real:: S_c_constant,S_p_constant

real:: dt                                           !time step
real:: time                                         !time
real:: error                                        !inertia error
real:: error_requirement                            !requirement of inertia error for convergence
logical:: convergence_index                         !true-convergence,false-not convergence

end module