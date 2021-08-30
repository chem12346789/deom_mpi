FROM continuumio/anaconda3:latest
RUN echo deb http://mirrors.ustc.edu.cn/debian buster main > /etc/apt/sources.list \
    && apt-get update\
    && apt-get install -y --allow-downgrades libssl-dev libssl1.1\
    && apt-get -y install tzdata\
    gcc make cmake g++ git curl libunwind-dev libblas-dev liblapack-dev liblapack3 libopenblas-base libopenblas-dev python python-numpy graphviz ghostscript nano gfortran autoconf\    
    libunwind8-dev \
    g++ wget gdb\
    cmake \
    libboost-all-dev \
    libevent-dev \
    libdouble-conversion-dev \
    libgoogle-glog-dev \
    libgflags-dev \
    libiberty-dev \
    liblz4-dev \
    liblzma-dev \
    libsnappy-dev \
    zlib1g-dev \
    binutils-dev \
    libjemalloc-dev \
    pkg-config \
    libunwind-dev \
    libelf-dev \
    libdwarf-dev
COPY gperftools /home/gperf
COPY json11 /home/json11
COPY jemalloc /home/jemalloc
COPY folly /home/folly
COPY googletest /home/googletest
COPY fmt /home/fmt
COPY eigen /home/eigen
RUN cd /home/gperf \
    && ./configure \
    && make -j$(nproc)\
    && make install \   
    && cd /home/json11 \
    && cmake . \
    && make -j$(nproc)\
    && make install \
    && cd /home/jemalloc \
    && . ./autogen.sh\
    && make -j$(nproc)\
    && make install\
    && rm -rf /home/gperf /home/json11 /home/jemalloc
RUN cd /home/googletest\
    && cmake .\
    && make -j6\
    && make install\
    && cd /home/fmt\
    && mkdir _build\
    && cd _build\
    && cmake ..\
    && make -j$(nproc)\
    && make install\
    && cd /home/folly\
    && mkdir _build && cd _build\
    && cmake configure -with-jemalloc ..\
    && make -j $(nproc)\
    && make install
RUN cd /home/eigen\
    && mkdir _build && cd _build\
    && cmake ..\
    && make install\
    && /sbin/ldconfig\
    && rm -rf /home/folly /home/googletest /home/fmt /home/eigen
