
# 使用基础镜像
ARG CUDA=12.2.0
FROM nvidia/cuda:${CUDA}-devel-ubuntu20.04

ENV DEBIAN_FRONTEND=noninteractive

# 安装必要的工具
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 复制打包的环境到容器
COPY inverse.tar.gz /opt/inverse.tar.gz

# 解压环境
RUN mkdir /opt/inverse && tar -xzf /opt/inverse.tar.gz -C /opt/inverse && \
    rm /opt/inverse.tar.gz

# 设置 PATH 环境变量并运行 conda-unpack
RUN export PATH=/opt/inverse/bin:$PATH && /opt/inverse/bin/conda-unpack

# 设置 PATH 环境变量
ENV PATH=/opt/inverse/bin:/opt/conda/bin:$PATH

# 设置工作目录
WORKDIR /workspace/esm-if

# 复制项目文件到容器
COPY ./esm-if1/ /workspace/esm-if

# 设置默认入口点，激活环境并运行脚本
#ENTRYPOINT ["/bin/bash", "-c", "source activate inverse && exec \"$@\""]

ENTRYPOINT ["/bin/bash", "-c", "source /opt/inverse/bin/activate inverse && /workspace/esm-if/run.sh \"$@\"", "--"]

#["/bin/bash", "-c", "source activate fastMSA && exec /workspace/DHR/fastMSA.sh \"$@\"", "--"]

# 定义容器启动时的命令
CMD ["AAAAAi", "1"]

