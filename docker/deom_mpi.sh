if [ -n "$1" ]; then
    docker run --pid=host --privileged -it --rm --workdir="/home/$1" --mount type=bind,source=$PWD/$1,target=/home/$1 deom_mpi:conda_dev
else
    echo "usage:deom <workdir>"
fi
