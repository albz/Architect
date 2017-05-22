# Guide for **background** namelist section <br/> background density tailoring

## Longitudinal profile Coordinates
To define the background *Longitudinal profile Coordinates* there are two equivalent options:
+ **z_coordinate_um**
+ **z_segment_length_um**

The *z_coordinate_um* corresponds to the laboratory global coordinate. It identifies each region where the background density is changing profile. The *z_segment_length_um* defines, instead, the length of each single segment. Note the z_segment_length_um is internally converted into z_coordinate_um, the code uses this latter quantity and its dimensionless form (*z_coordinate*).

## Transverse Density profile

In order to set the radial profile of the background plasma it is possible to choose between some schemes by setting:

+ **order\_radial=0**
+ **order\_radial=1**
+ **order\_radial=2**
+ **order\_radial=3**
+ **order\_radial=4**
+ **order\_radial=5**

In the following the radial coordinate is indicated with the letter $r$.

### order\_radial=0
Uniform density in $r$

### order\_radial=3
The profile starts from a fixed value at r=0 and decays with a $cos^2$ shape reaching zero  at $r=$**radius\_um**.

### order\_radial=4
The profile has a cos<sup>2</sup> shape shrunk vertically by a factor **perturbation\_amplitude** and shifted upwards by a **certain value**. The profile drops to 0 with a step at $r=$**radius\_um**.

### order\_radial=5
The profile has a flat-top with a **certain hight** from $r=0$ to $r=$**radius\_internal\_um** then, from **radius\_internal\_um** to **radius\_um** the profile decays with a $cos^2$ shape.

## Longitudinal Density profile
In order to set the radial profile of the background plasma it is possible to choose between some schemes by setting:
+ **bck_plasma%order_logitudinal=0**
+ **bck_plasma%order_logitudinal=1**
+ **bck_plasma%order_logitudinal=2 (and -2)**
+ **bck_plasma%order_logitudinal=3**

In the following the longitudinal coordinate is indicated with the letter $z$.

### order_logitudinal=0
Uniform profile in $z$. It is necessary to set the value of the (relative) plasma density at both the boundaries of the subdomain where the flat-top profile is desired. If $k$ is the subdomain to be set, bck_plasma%n_over_n0(k) and bck_plasma%n_over_n0(k+1) must be both equal to the wanted relative density.

### order_logitudinal=1
The profile is made by a liner ramp between the plasma density values which are set at the boundary of the subdomain. If $k$ is the subdomain to be set, bck_plasma%n_over_n0(k) and bck_plasma%n_over_n0(k+1) must be equal to the wanted relative density at the borders of the linear ramp.

### order_logitudinal=2
Concave downward parabolic profile. The value of the plasma density set at boundaries of the subdomain (by assigning a value to bck_plasma%n_over_n0(k) and bck_plasma%n_over_n0(k+1), for the $k$-th longitudinal subdomain) are connected with a segment of parabola which is always concave downward and its vertex is on the boundary with higher density value.

### order_logitudinal=-2
Convex upward parabolic profile. The value of the plasma density set at boundaries of the subdomain (by assigning a value to bck_plasma%n_over_n0(k) and bck_plasma%n_over_n0(k+1), for the $k$-th longitudinal subdomain) are connected with a segment of parabola which is always convex upward and its vertex is on the boundary with lower density value.

### order_logitudinal=3
The profile has a $cos^2$ shape shrunk and shifted in order to connect the value of plasma density set at the subdomain borders.
