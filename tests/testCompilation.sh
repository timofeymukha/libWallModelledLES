if [ $1 == "4.1" ]
then
    docker start of_4
    docker exec -w $PWD/.. of_4 /bin/bash -c ". /opt/openfoam4/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-4.1/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake"
fi
if [ $1 == "5.0" ]
then
    docker start of_5
    docker exec -w $PWD/.. of_5 /bin/bash -c ". /opt/openfoam5/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-5/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake"
fi
if [ $1 == "6.0" ]
then
    docker start of_6
    docker exec -w $PWD/.. of_6 /bin/bash -c ". /opt/openfoam6/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-6/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake"
fi
if [ $1 == "7.0" ]
then
    docker start of_7
    docker exec -w $PWD/.. of_7 /bin/bash -c ". /opt/openfoam7/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-7/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake"
fi
if [ $1 == "8.0" ]
then
    docker start of_8
    docker exec -w $PWD/.. of_8 /bin/bash -c ". /opt/openfoam8/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-7/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake"
fi
if [ $1 == "9.0" ]
then
    docker start of_9
    docker exec -w $PWD/.. of_9 /bin/bash -c ". /opt/openfoam9/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-7/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v30plus" ]
then
    docker start of_v30plus
    docker exec -w $PWD/.. of_v30plus /bin/bash -c "source /opt/OpenFOAM/OpenFOAM-v3.0+/etc/bashrc && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v1606" ]
then
    docker start of_v1606
    docker exec -w $PWD/.. of_v1606 /bin/bash -c "source /opt/OpenFOAM/setImage_v1606+ && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v1612" ]
then
    docker start of_v1612
    docker exec -w $PWD/.. of_v1612 /bin/bash -c "source /opt/OpenFOAM/setImage_v1612+ && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v1706" ]
then
    docker start of_v1706
    docker exec -w $PWD/.. of_v1706 /bin/bash -c "source /opt/OpenFOAM/setImage_v1706.sh && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v1712" ]
then
    docker start of_v1712
    docker exec -w $PWD/.. of_v1712 /bin/bash -c "source /opt/OpenFOAM/setImage_v1712.sh && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v1806" ]
then
    docker start of_v1806
    docker exec -w $PWD/.. of_v1806 /bin/bash -c "source /opt/OpenFOAM/setImage_v1806.sh && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v1812" ]
then
    docker start of_v1812
    docker exec -w $PWD/.. of_v1812 /bin/bash -c "source /opt/OpenFOAM/setImage_v1812.sh && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v1906" ]
then
    docker start of_v1906
    docker exec -w $PWD/.. of_v1906 /bin/bash -c "source /opt/OpenFOAM/setImage_v1906.sh && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v1912" ]
then
    docker start of_v1912
    docker exec -w $PWD/.. of_v1912 /bin/bash -c "source /opt/OpenFOAM/setImage_v1912.sh && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v2006" ]
then
    docker start of_v2006
    docker exec -w $PWD/.. of_v2006 /bin/bash -c "source /opt/OpenFOAM/setImage_v2006.sh && ./Allwclean && ./Allwmake"
fi
if [ $1 == "v2012" ]
then
    docker start of_v2012
    docker exec -w $PWD/.. of_v2012 /bin/bash -c "source /opt/OpenFOAM/setImage_v2012.sh && ./Allwclean && ./Allwmake"
fi
if [ $1 == "all" ]
then
    mkdir logs
    docker start of_4
    docker exec -w $PWD/.. of_4 /bin/bash -c ". /opt/openfoam4/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-4.1/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake &> tests/logs/4.log"
    docker start of_5
    docker exec -w $PWD/.. of_5 /bin/bash -c ". /opt/openfoam5/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-5/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake &> tests/logs/5.log"
    docker start of_6
    docker exec -w $PWD/.. of_6 /bin/bash -c ". /opt/openfoam6/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-6/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake &> tests/logs/6.log"
    docker start of_7
    docker exec -w $PWD/.. of_7 /bin/bash -c ". /opt/openfoam7/etc/bashrc && export FOAM_USER_LIBBIN=/home/timofey/OpenFOAM/timofey-6/platforms/linux64GccDPInt32Opt/lib && ./Allwclean && ./Allwmake &> tests/logs/7.log"
    docker start of_v30plus
    docker exec -w $PWD/.. of_v30plus /bin/bash -c "source /opt/OpenFOAM/OpenFOAM-v3.0+/etc/bashrc && ./Allwclean && ./Allwmake &> tests/logs/3+.log"
    docker start of_v1606
    docker exec -w $PWD/.. of_v1606 /bin/bash -c "source /opt/OpenFOAM/setImage_v1606+ && ./Allwclean && ./Allwmake &> tests/logs/1606.log"
    docker start of_v1612
    docker exec -w $PWD/.. of_v1612 /bin/bash -c "source /opt/OpenFOAM/setImage_v1612+ && ./Allwclean && ./Allwmake &> tests/logs/1612.log"
    docker start of_v1706
    docker exec -w $PWD/.. of_v1706 /bin/bash -c "source /opt/OpenFOAM/setImage_v1706.sh && ./Allwclean && ./Allwmake &> tests/logs/1706.log"
    docker start of_v1712
    docker exec -w $PWD/.. of_v1712 /bin/bash -c "source /opt/OpenFOAM/setImage_v1712.sh && ./Allwclean && ./Allwmake &> tests/logs/1712.log"
    docker start of_v1806
    docker exec -w $PWD/.. of_v1806 /bin/bash -c "source /opt/OpenFOAM/setImage_v1806.sh && ./Allwclean && ./Allwmake &> tests/logs/1806.log"
    docker start of_v1812
    docker exec -w $PWD/.. of_v1812 /bin/bash -c "source /opt/OpenFOAM/setImage_v1812.sh && ./Allwclean && ./Allwmake &> tests/logs/1812.log"
    docker start of_v1906
    docker exec -w $PWD/.. of_v1906 /bin/bash -c "source /opt/OpenFOAM/setImage_v1906.sh && ./Allwclean && ./Allwmake &> tetts/logs/1906.log"
fi
