import numpy as np
from nemotoc.py_transform.tom_sum_rotation import tom_sum_rotation
from nemotoc.py_transform.tom_pointrotate import tom_pointrotate
from nemotoc.py_transform.tom_angular_distance import tom_angular_distance

def tom_calcPairTransForm(pos1, ang1, pos2, ang2, dMetric = 'exact',symmetry=1):
    '''
    TOM_CALCPAIRTRANSFORM calculates relative transformation between two poses
                 
    [posTr,angTr]=tom_calcPairTransForm(pos1,ang1,pos2,ang2)
    
    PARAMETERS
    
    INPUT
        pos1                 position 1
        ang1                 angle1 in Z-X-Z
        pos2                 posion2
        ang2                 angle2 in Z-X-Z
        dMetric              ('exact') at the moment only exact implemented
        symmetry             the rotational symmetry
    
    OUTPUT
        posTr               transformation vector between two points (the coordinates(x-y-z) of one pose coordinate system) 
                            one-dimensition array
                            
        angTr               transformation angle between two points
                            one-dimensition array
        
        lenPosTr            length of transformation vector     
        lenAngTr            angular distance from [0 0 0] to angTr (using quaternions)
    
    
    EXAMPLE
       [pos,rot]=tom_calcPairTransForm(np.array([1, 1, 1]),np.array([0, 0, 10]),np.array([2, 2, 2]),np.array([0, 0, 30]));
      
    
    REFERENCES    
    '''
    if dMetric == 'exact':
        angRot = 360/float(symmetry)
        ang1Inv = np.array([-ang1[1],-ang1[0],-ang1[2]])
        ang2Inv = np.array([-ang2[1],-ang2[0],-ang2[2]])
        angTrs = [] # the angle transforms of symmetry. 
        posTrs = []
        lenPosTrs = []
        lenAngTrs = []
        # add M0 rotation
        for i in range(symmetry):
            compare_array = np.zeros([3,3])
            compare_array[0,:] = ang2
            compare_array[1,:] = np.array([0,-angRot*i,0])
            compare_array[2,:] = ang1Inv
            #calculate euler angles of relative rotation
            angTr, _, _ = tom_sum_rotation(compare_array, np.zeros([3,3]))
            angTrs.append(angTr)
            #calclulate relative coordinates of pos2 in the ang1-po1 coordinate system
            pos2Rel = pos2-pos1
            tmp_ang,_,_ = tom_sum_rotation(compare_array[1:,:],np.zeros([2,3]))
            posTr = tom_pointrotate(pos2Rel, tmp_ang[0], tmp_ang[1], tmp_ang[2])
            posTrs.append(posTr)
            #calculte eluer distance of posTr
            lenPosTr = np.linalg.norm(posTr)
            lenPosTrs.append(lenPosTr)
            #calculate the quaternions, from zero angle to relative rotation angle
            lenAngTr = tom_angular_distance(np.zeros(3),angTr)
            lenAngTrs.append(lenAngTr)
        index1 = lenPosTrs.index(min(lenPosTrs))
        index2 = lenAngTrs.index(min(lenAngTrs))
    return posTrs[index1], angTrs[index2], lenPosTrs[index1], lenAngTrs[index2]   