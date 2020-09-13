# SR_ROBO

<p float="center">
  <img src="/IMAGES/uni_rrobot_u10.png" width="100%%" title=" "/>
</p>

MATLAB Robotics custom toolbox and examples, linear algebra functions, manipulator's direct kinematic tools, inverse direct kinematic, numerical and analytical methods, differential kinematic, inverse differential kinematic, trajectory planning functions at cartesian and joint space levels, some controle scheme and some basic feedforward-feedback control loops. 
Here is a basic example on what you can achieve with this toolbox:

**Anthropomorphic arm with spherical wrist (6R) performing a circular path with end effector changing orientation over quintic time law profile.**

**Manipulator:**
<p  float="center">
  <img src="/IMAGES/6R.png" width="50%%" title=" "/>
</p>

    joint_i     alpha_i      a_i       d_i      theta_i
    _______    _________    ______    ______    _______

       1       {'pi/2' }    {'0' }    {'0' }    {'q1'} 
       2       {'0'    }    {'a2'}    {'0' }    {'q2'} 
       3       {'pi/2' }    {'a3'}    {'0' }    {'q3'} 
       4       {'-pi/2'}    {'0' }    {'d4'}    {'q4'} 
       5       {'pi/2' }    {'0' }    {'0' }    {'q5'} 
       6       {'0'    }    {'0' }    {'0' }    {'q6'} 
       
  **Trajectory:**
 <p float="center">
  <img src="/IMAGES/traj.png" width="50%%" title=" "/>
</p>
