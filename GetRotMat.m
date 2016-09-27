function RT_abc=GetRotMat(xztheta, loop_t, loop_phi, loop_v)

new_states=C_Predict_P(xztheta,loop_t,loop_phi,loop_v);
d_theta_hump=-new_states(3)+pi/2;
R_hump=[cos(d_theta_hump) -sin(d_theta_hump); sin(d_theta_hump) cos(d_theta_hump)];
T_hump=-R_hump*[new_states(1); new_states(2)];
RT_abc=[R_hump T_hump; 0 0 1];

end