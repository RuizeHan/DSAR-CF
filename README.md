# DSAR-CF
Public code of DSAR-CF (Dynamic Saliency-Aware Regularization for Correlation Filter based Object Tracking), published in TIP 2019.
```
@article{han2019dynamic,
  title={Dynamic Saliency-Aware Regularization for Correlation Filter based Object Tracking}, 
  author={Feng, Wei and Han, Ruize and Guo, Qing and Zhu, Jianke and Wang, Song},  
  year={2019},  
  booktitle={IEEE Transactions on Image Processing}
}
```

## Method

![example](https://github.com/HanRuize/DSAR-CF/blob/master/figs/example.png)

In the proposed DSAR-CF, we first introduce object saliency information into the regularization weight map to highlight the appearance of the player as well as suppressing the background information around the player at the first frame.  
We then propose a strategy to dynamically update the regularization weight map to reflect the player’s shape variations in the subsequent frames.  
We then develop a level-set algorithm to iteratively optimize the regularization weight map in each frame. 

## Experiments
The experimental results on OTB-2013/2015:  

![res](https://github.com/HanRuize/DSAR-CF/blob/master/figs/res.png)
