## some (work in progress) background density tailoring 

In order to set the radial profile of the background plasma it is possible to choose between some schemes by setting:

+ **order\_radial::0**
+ **order\_radial::1**
+ **order\_radial::2**
+ **order\_radial::3**
+ **order\_radial::4**
+ **order\_radial::5**

In the following the radial coordinate is indicated with the letter _r_.

### order\_radial::0
Uniform density in _r_

### order\_radial::3
The profile starts from a fixed value at r=0 and decays with a cos<sup>2</sup> shape reaching zero  at _r_=**radius\_um**.

### order\_radial::4
The profile has a cos<sup>2</sup> shape shrunk vertically by a factor **perturbation\_amplitude** and shifted upwards by a **certain value**. The profile drops to 0 with a step at _r_=**radius\_um**.

### order\_radial::5
The profile has a flat-top with a **certain hight** from _r_=0 to r=**radius\_internal\_um** then, from **radius\_internal\_um** to **radius\_um** the profile decays with a cos<sup>2</sup> shape.
