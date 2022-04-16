import numpy as np
# u v pt ps delta_u delta_v delta_pt delta_ps

def Colorrange(testcase):
    if testcase == 'u':
        colorrange = np.array([[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10],[-20,20], [-20,20]])
    elif testcase == 'v':
        colorrange = np.array([[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10],[-20,20], [-20,20]])
    elif testcase == 't':
        colorrange = np.array([[130, 140],[150,170], [190, 210], [190, 210], [210, 250], [270, 285], [280,290]])
    elif testcase == 'pt':
        colorrange = np.array([[280, 300],[280,300], [300, 350], [300, 350], [350, 400], [350, 400], [440,500]])
    elif testcase == 'ps':
        colorrange = np.array([[95000, 105000],[95000, 105000],[95000, 105000],[95000, 105000],[95000, 105000],[95000, 105000],[95000, 105000]])
    return colorrange
           
def Deltacolorrange(testcase):
    if testcase == 'u':
        deltacolorrange = np.array([[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10],[-20,20], [-20,20]])
    elif testcase == 'v':
        deltacolorrange = np.array([[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10],[-20,20], [-20,20]])
    elif testcase == 't':
        deltacolorrange = np.array([[130, 140],[150,170], [190, 210], [0, 25], [210, 250], [270, 285], [280,290]])
    elif testcase == 'pt':
        deltacolorrange = np.array([[130, 140],[150,170], [190, 210], [0, 25], [210, 250], [270, 285], [280,290]])
    elif testcase == 'ps':
        deltacolorrange = np.array([[-100, 100],[95000, 105000],[95000, 105000],[95000, 105000],[95000, 105000],[95000, 105000],[95000, 105000]])
    return deltacolorrange
    
     #设定变量的range(每4层一个range)：最小值和最大值
