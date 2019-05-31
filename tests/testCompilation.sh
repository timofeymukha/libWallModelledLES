if [ $1 == "v30plus" ]
then
    docker start of_v30plus
    docker exec -w $PWD/.. of_30plus /bin/bash -c "source /opt/OpenFOAM/setImage_v1706.sh && ./Allwclean && ./Allwmake"
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

#docker exec -w $PWD/.. of_v1706 /bin/bash -c "source /opt/OpenFOAM/setImage_v1706.sh && ./Allwclean && ./Allwmake"
#docker exec -w $PWD/.. of_v1712 /bin/bash -c "source /opt/OpenFOAM/setImage_v1712.sh && ./Allwclean && ./Allwmake"
#docker exec -w $PWD/.. of_v1806 /bin/bash -c "source /opt/OpenFOAM/setImage_v1806.sh && ./Allwclean && ./Allwmake"
#docker exec -w $PWD/.. of_v1812 /bin/bash -c "source /opt/OpenFOAM/setImage_v1812.sh && ./Allwclean && ./Allwmake"

