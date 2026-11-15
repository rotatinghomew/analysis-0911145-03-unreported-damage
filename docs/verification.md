# 运行与验证

现有业务预览、资料校验及真实 HTTP 传输由 `python -m unittest discover -s tests -v` 验证；`python -m app.demo` 临时启动服务，通过消费端读取旧业务预览、新资料描述及 `POST /pooled-quality` 的真实分析输出。固定依赖的安装和容器准备方法见 README。

`RUN_ACCEPTANCE=1 python -m unittest tests.test_acceptance -v` 验证分析消费接口的基本字段与信息不足契约（接口已接通，默认仍由环境变量控制启用）。该测试不是算法正确性的证明；算法正确性由 `tests/test_analysis.py` 覆盖：确定性仪器与点校准两个可手算小例（联合似然、共享果实、计数配置后验与手算一致）、重复入管与追加管的设计核算、缺校准管剔除且不借用其他大小、中点区间落于角点包络内的敏感性回归、复检额度与唯一性契约，以及真实 HTTP 传输。演示按资料名请求真实分析路由，不返回静态结论；正式展示不安装依赖，资料规模在数秒内结束。
