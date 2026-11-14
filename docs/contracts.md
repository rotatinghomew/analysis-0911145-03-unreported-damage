# 产品数据与调用约定

population 是发货总体，不含真实货损隐藏字段。lists 中来源为 service、carrier、callback。callback 的 outreach 包含抽样概率、是否接通，没接通不等于无货损；主动渠道无已知纳入概率。confirmed_links 是可直接合并的身份，candidate_links 仅表示可能相同，可同时包含互斥候选；不能把候选分值当真实概率。POST /coverage {case, stratum?} 返回 {snapshot, observed, estimates, assumptions, diagnostics}。observed 报告已确认可见事件；estimates 须区分点、区间、敏感性情景及不可识别。不得把总体未捕获者写成确认为无货损。样本皆为本地生成的业务记录。

当前 GET /cases、POST /preview 和 GET /health 已可用；client.analysis 是待接入的业务消费调用，新分析路由尚不存在。原始数据通过 app.store 读取，snapshot 为输入字节的 SHA-256。新增接口须保留现有预览能力和 JSON 有限数值约定。

observed.confirmed_events 是至少被一条确认名单捕获的不同订单数，而非记录条数。repeated 含同一渠道的重复记录，应去重后参与重叠。estimates 可按不同依赖假设给区间，但必须标明纳入的身份情景、分层总体及未回答回访的处理。
