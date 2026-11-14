# 产品数据与调用约定

specimens 的 id 是同一果实，lot 为采购批次，inclusion_probability 是最初从 lot_population 抽中该果的概率。pools 的 members 引用果实，positive 是混样仪器二元读数而非逐果标签；stage/selected_after 表示根据先前管结果追加检测的策略。calibration 按 assay、pool_size、bad_count 给出阳性概率区间；每一行区间不是额外随机独立的仪器来源。未提供某大小的校准不得静默沿用其他大小。新接口 POST /pooled-quality 接收 {case, audit_budget?}，返回 {snapshot, lots, identification, followup, diagnostics}。lots 区分样本缺陷数和批次比例，范围须注明可信或可行含义；identification 解释仍相容的配置，followup 受额度及实际可复检果实约束。

现有 GET /health、GET /cases 和 POST /preview 保留原有资料预览。新增 GET /study-cases 和 POST /study-preview {case} 已可用，只校验与描述 data/study 的资料；推断消费路由尚未实现。app.client.analysis 指向上文约定的新消费接口。所有 snapshot 都是相应原始资料字节的 SHA-256，输出仅使用有限 JSON 数值，无法判断用明确状态与 null 表达，不得返回 NaN。资料只供本地分析，不连接外部服务。

lots 是记录列表，每项含 lot、sample_defects 和 prevalence；sample_defects 为被抽中果实的缺陷数范围或 null，prevalence 为批次比例范围或 null。followup 含 selected_specimens 列表，不得重复消耗同一果实的额度。
