# DEOM

## FILE  

This includes two-part (bose and fermi) of [HEOM](https://en.wikipedia.org/wiki/Hierarchical_equations_of_motion) (Hierarchical equations of motion):  

1. bose: boson bath and boson system
2. fermi: fermi bath and fermi system
3. docker: docker image builder
   1. dev:  development environment with python, it will cause image size very huge.
   2. deom_mpi.sh: a bash shell help you to run the application.

and as the name of the folder indicate, 1d-corr means in that folder, you can calculate the equilibrium state of a given system, and calculate the correlation function of that. The sto_quad folder contains a stochastic algorithm to calculate quad system-bath coupling.

## INSTALL

You should consider trying [docker](https://www.docker.com/), docker image builder is in /docker folder. And note if you are in USTC, a harbor service that contains this image can be provided, you can contact us.  
In the /docker folder, we only provide Dockerfile, you have to download those software and modify version number.

1. [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
2. [folly](https://github.com/facebook/folly/tree/master/folly)
   1. [googletest-release](https://github.com/google/googletest/releases)
   2. [fmt](https://github.com/fmtlib/fmt)
   3. [gperftools](https://github.com/gperftools/gperftools)
   4. [jemalloc](https://github.com/jemalloc/jemalloc)(get into trouble? [click it!](https://github.com/facebook/folly/issues/943))
3. [json11](https://github.com/dropbox/json11)

DO NOT forget to modify the version number! Then you can build the image with `docker build  -t deom_mpi:conda_dev .`, and you can run a container by deom_mpi.sh in /docker directory. You can change the deom_mpi:conda_dev to IMAGE: TAG you like but DO NOT forget to modify deom_mpi.sh after doing that. In docker container just build the application you need and then run it!

## BRIEF

This is a HEOM calculator that utilizes the hashmap to do parallelize filter. Technical details is very simple:
There are three things we need to do:

1. Store all density matrices
2. Look for the density matrix used
3. Calculate the super operator and calculate the $\dot rho^(n)_{/bf n}$

in step 1 and step 2 we will use an atomic hash map to do search things. In step 3, we can do it by using an embarrassingly parallel algorithm.

## NOTE

We WON NOT provide input file samples in some applications due to that part of the work is not complete.

## translation

English document is up to here, 下面是中文文档。

## 文件

这包括 [HEOM](https://en.wikipedia.org/wiki/Hierarchical_equations_of_motion)（级联运动方程）的两部分（玻色和费米）：

1.bose：玻色子浴和玻色子系统
2. fermi：费米浴和费米系统
3. docker: docker 镜像构建器

   1. dev：使用python开发环境，会导致图片体积很大。
   2. deom_mpi.sh：一个帮助你运行应用程序的bash shell。

正如文件夹名称所示，1d-corr 表示在该文件夹中，您可以计算给定系统的平衡状态，并计算其相关函数。 sto_quad 文件夹包含一个随机算法来计算二次系统-环境耦合。

## 安装

您应该考虑尝试 [docker](https://www.docker.com/)，docker image builder 位于 /docker 文件夹中。请注意，如果您在中国科学技术大学，可以提供包含此图像的harbor服务，您可以联系我们。
在 /docker 文件夹中，我们只提供 Dockerfile，您需要下载那些软件并修改版本号。

1. [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
2. [folly](https://github.com/facebook/folly/tree/master/folly)
   1. [googletest-release](https://github.com/google/googletest/releases)
   2. [fmt](https://github.com/fmtlib/fmt)
   3. [gperftools](https://github.com/gperftools/gperftools)
   4. [jemalloc](https://github.com/jemalloc/jemalloc)(遇到麻烦？[点击它！](https://github.com/facebook/folly/issues/943))
3. [json11](https://github.com/dropbox/json11)

不要忘记修改版本号！然后你可以使用`docker build -t deom_mpi:conda_dev .`构建镜像，你可以通过/docker目录下的deom_mpi.sh运行一个容器。您可以将 deom_mpi:conda_dev 更改为您喜欢的IMAGE:TAG，但不要忘记在执行此操作后修改 deom_mpi.sh。在 docker 容器中，只需构建您需要的应​​用程序，然后运行它！

## 简介

这是一个 HEOM 计算器，它利用 hashmap 进行并行化过滤。技术细节很简单：
我们需要做三件事：

1. 存储所有密度矩阵
2. 寻找使用的密度矩阵
3. 计算超算符，计算$\dot rho^(n)_{/bf n}$

在第1步和第2步中，我们将使用原子哈希映射来进行搜索。在第3步中，我们可以使用令人尴尬的并行算法来完成。

## 注意

由于部分工作未完成，我们不会在某些应用程序中提供输入文件示例。
