# SR_ROBO

**AUTHOR:** 

*Simone Rossetti, Rome, August 2020*

[Istitutional](mailto:rossetti.1900592@studenti.uniroma1.it)

[Private](mailto:simone.rossetti@live.com)

**NB: The use of this material is free as long as the content is not altered and the name of the author and the link to this repo are always cited.**

<p float="center">
  <img src="/IMAGES/uni_rrobot_u10.png" width="70%%" title=" "/>
</p>

MATLAB Robotic custom toolbox and examples. 
It contains: linear algebra functions, manipulator's direct kinematic tools, inverse direct kinematic, numerical and analytical methods, differential kinematic, inverse differential kinematic, trajectory planning functions at cartesian and joint space levels, some controle scheme and some basic feedforward-feedback control loops. 
Here is a basic example on what you can achieve with this toolbox:

**Anthropomorphic arm with spherical wrist (6R) performing a circular path with end effector changing orientation over quintic time law profile.**

**DH table:**

    joint_i     alpha_i      a_i       d_i      theta_i
    _______    _________    ______    ______    _______

       1       {'pi/2' }    {'0' }    {'0' }    {'q1'} 
       2       {'0'    }    {'a2'}    {'0' }    {'q2'} 
       3       {'pi/2' }    {'a3'}    {'0' }    {'q3'} 
       4       {'-pi/2'}    {'0' }    {'d4'}    {'q4'} 
       5       {'pi/2' }    {'0' }    {'0' }    {'q5'} 
       6       {'0'    }    {'0' }    {'0' }    {'q6'} 
       
  **Manipulator:**
<p  float="center">
  <img src="/IMAGES/6R.png" width="50%%" title=" "/>
</p>

**Time law:**
 <p float="center">
  <img src="/IMAGES/time.png" width="50%%" title=" "/>
</p>

  **Trajectory:**
 <p float="center">
  <img src="/IMAGES/traj.png" width="50%%" title=" "/>
</p>

**End effector velocity:**
 <p float="center">
  <img src="/IMAGES/vel.png" width="50%%" title=" "/>
</p>

**End effector rotation (axis/angle):**
 <p float="center">
  <img src="/IMAGES/ee_angle.png" width="50%%" title=" "/>
</p>

**End effector angular velocity:**
 <p float="center">
  <img src="/IMAGES/ee_vel.png" width="50%%" title=" "/>
</p>

**Result:**
 <p float="center">
  <img src="/IMAGES/res.png" width="45%%" title=" "/>
    <img src="/IMAGES/res2.png" width="45%%" title=" "/>
</p>
