1.基于 qcat库，该库在github中管理，而pland也使用。最初是基于同样的原版本修改而来，后期github可能会进一步更改。
并以GitHub为准而pland跟进，尽量不改动。作为子模块使用



2.基于 分配算法对比，在pland中有使用，最初是在航空学报录用文件夹，改函数名后到pland/alloc_compa使用。
保留航空学报的版本将不再改动，而基于pland的版本重新成立新文件夹进行git管理。pland将基于此github版。主要是ac.slx，全部复制了11个文件


3.基于 allocation m2C，使用几何法观察可达集。实际上只用一个test和constraint函数和twin.slx。


宗旨：用独立的仓库管理(实际上都从当前的pland的ASUS_d分支， 8a49c8570c6e1ed4ea218c7c44a3319d2d84a537开始)，
pland跟进，旧文件（TDF、航空学报论文）的内容不改动。
