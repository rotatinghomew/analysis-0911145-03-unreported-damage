# 重叠混样与稀释误差下的果品批次缺陷比例推断

本地分析工作台保留已有业务资料与预览能力，并接入了混样筛查的联合推断。`data/study` 保存本地合成资料，`app/study.py` 解释资料结构，`app/pooled.py` 实现推断：同一果实的缺陷状态联合进入其所有混管的似然（不按管独立伯努利、重复入管不计新增样本量），校准区间按系统性误差在三个角点（下界/中点/上界）上精确计算并报告包络，缺校准的混管读数被剔除并列入诊断，追加拆样作为联合似然因子而非新随机样本。描述性预览（`POST /study-preview`）只描述观测，推断结果以 `POST /pooled-quality` 为准。

数据均为本地合成，不含真实姓名、电话号码或生产系统凭据。`standard` 是常规资料，`edge` 用于检查信息不足或约束冲突。它们只描述观测，不预设分析结论。

## 本地使用

Python 3.11，依赖版本固定在 `requirements.txt`。NumPy/SciPy 是分析工作台的数值工具；预览逻辑主要使用标准库。准备阶段安装依赖，仓库不在启动或测试时下载任何资源。

支持 Python 3.11 系列的补丁版本；Docker 基础镜像通过 SHA256 固定为 Python 3.11.14，依赖清单保持精确版本。相同依赖还需在非 root 的 Python 3.11.2 Linux 环境中验证，应用不依赖特定补丁版本。

```sh
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python -m app.demo
python -m app.server
```

服务监听 `127.0.0.1:8080`，也可通过 `PORT` 调整。预览请求由 `app/client.py` 发起；`app.demo` 在临时端口启动真实服务并完整打印现有业务输出，不是录制日志的回放。服务接口和数据文件的可见性约定不能互相替代。

## 可复现环境

```sh
docker build -t local/unreported-damage:baseline .
docker run --rm --network none local/unreported-damage:baseline python -m unittest discover -s tests -v
docker run --rm --network none local/unreported-damage:baseline
```

准备阶段允许构建镜像和预装依赖；正式演示不进行构建或安装。演示通过实际 HTTP 消费调用展示输入口径、分析结果和限制，保留完整输出；小样本在数秒内结束。功能测试启用方式：`RUN_ACCEPTANCE=1 python -m unittest tests.test_acceptance -v`；分析正确性回归（小规模可手算核对、敏感性、信息截止）见 `tests/test_analysis.py`，默认随测试套件运行。
