# 行为验证

- 多来源捕获重叠分析
- 依赖与抽访设计下的漏报识别
- 身份歧义的估计传播
- 分层诊断与新增调查设计

tests/test_domain.py 覆盖现有业务口径；tests/test_transport.py 通过真实本地 HTTP 验证数据到消费端。tests/test_acceptance.py 是新接口的公开契约测试，完成新能力后以 RUN_ACCEPTANCE=1 启用。它不替代统计正确性、不可识别性、组合约束和敏感性回归；请按 contracts.md 补足针对算法设计的测试。新接口不得通过固定 fixture 名称返回预制结论。
