运行：./run.sh




run.sh：运行脚本
for test：仅针对已存入的真实场层[1,32,40,59]测试用
for use：针对集成测试用，配合consumer脚本中set space部分注释掉的配置

consumer-mpi-multi-real.py：主脚本
start_change->end_change：可修改部分

testplotmulti.py：绘图脚本

cal_indicator.py：指标计算

colorrangemulti.py：色标给定
此处暂将60层定为每10层一个色标,可根据集成测试更改

