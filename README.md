# stereovision-on-road-detection

This repository is open to public, any one wants to implement road lane line detection, can freely use it. 
It integrates a new lane line detection algorithm with other lane marking detectors to identify the correct lane line markings. It also fits multiple road models to improve accuracy. An effective stereo 3D reconstruction method is proposed to estimate vehicle localization. The estimation consistency is further guaranteed by a new particle filter framework, which takes vehicle dynamics into account. Experiment results based on image sequences taken under different visual conditions showed that the proposed system can identify the lane line markings with 98.6% accuracy. The maximum estimation error of the vehicle distance to lane lines is 16cm in daytime 26cm at night, and the maximum estimation error of its moving direction respected to road tangent is 0.06rad in daytime and 0.12rad at night. 
For more details and for citation, please refer to the following:

@article{du2016comprehensive,
  title={Comprehensive and practical vision system for self-driving vehicle lane-level localization},
  author={Du, Xinxin and Tan, Kok Kiong},
  journal={IEEE transactions on image processing},
  volume={25},
  number={5},
  pages={2075--2088},
  year={2016},
  publisher={IEEE}
}
