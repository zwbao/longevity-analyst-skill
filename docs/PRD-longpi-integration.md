# PRD：longpi（dsh 插件）× longevity-analyst 适配（草案，待确认）

日期：2026-09-30。依据：本机部署 dsh 0.1.5 + longpi 0.6.3 + Mirobody 1.5.2（原生）+ longevity-skills 2026.40.0 +
longevity-analyst，使用 DeepSeek（deepseek-v4-pro）和虚拟会员「李明华」做的实测。

## 1. 最终形态

会员只和 dsh 上的 longpi 页面打交道。longpi 负责日常：数据录入（Mirobody）、方案、打卡、随访、提醒。
longevity-analyst 负责深度分析：多组学、器官体检表、洞见层（基因对照、人群位置、因果推算、问题看板）、报告和孪生快照。
两者共用同一份 longevity-skills 方法库，由 longpi 安装器放在 `~/longpi/longevity-skills`。

| 环节 | 现在 | 目标 |
|---|---|---|
| 数据进来 | analyst 只读原始文件夹；longpi 只读 Mirobody | analyst 能直接从 Mirobody 拿会员已确认的化验、可穿戴日汇总和上传的报告文件 |
| 触发 | 会员在对话里点名 longevity-analyst | longpi 首页或「健康」页有「深度分析」入口，一键发起，进度可见 |
| 结果回去 | analyst 只产出 deliver/ 下的 report 和 twin | 结果写回 longpi：报告与问题看板进「健康」页新 tab；读数进「指标」；干预方案转成 longpi 方案草稿，由会员读回确认；复测计划进 longpi 随访 |
| 复测 | analyst 手动 twin compare | longpi 到期提醒 → 新一轮分析 → 自动与上期孪生快照对比 |

## 2. 实测发现的问题（按优先级）

| 级别 | 问题 | 证据 | 归属 |
|---|---|---|---|
| P0 | dsh 默认把会话日志随请求上传给 DeepSeek（session-log-deepseek），健康和基因数据会进 DeepSeek 日志 | dsh 文档 python-sdk.md；本机已在 home patch 关闭 | longpi 安装器默认写入关闭 |
| P0 | longpi 安装器找 LOINC 词表的路径已过时（Mirobody 把它移到 `res/loinc/`），导致 Mirobody 部署失败 | install.log；手动指定 LONGPI_LOINC_URL 后通过 | longpi |
| P1 | 国内网络下 Mirobody Docker 镜像构建失败（apt 源经代理不稳定），`thetahealth/mirobody:1.5.3` 未发布到 Docker Hub | 两次构建失败日志 | longpi 安装器增加原生运行模式（已验证：pg 容器 + 本机 `mirobody serve/worker`） |
| P1 | analyst 与 longpi/Mirobody 零接口 | 代码 grep 两边互不引用 | 本 PRD 第 3 节 |
| P1 | 方法库里 `testis-transcriptomic-atlas-lifespan` 被标为 verified 而自动运行，输出「年龄段 50 多岁」这种无关结果 | 李明华实测 readouts | longevity-skills |
| P2 | longpi 文档写 14 个工具，代码注册约 41 个 | 调研报告 | longpi |
| P2 | Mirobody 本地无邮件服务时只有预置演示邮箱能登录，新会员需 `/password/register` | 本机实测 | 文档 |

## 3. 适配改造

### 3.1 analyst 侧（本仓库）
1. `la.py intake --from-mirobody <mcp-url>`：通过会员个人 MCP（只读）导出已确认化验（带 LOINC、单位、参考区间、采样日期）、
   可穿戴日汇总和上传原件到 workspace 的 raw 目录，再走原有 intake；来源标记为 mirobody，行确认沿用 `labs confirm`。
2. `la.py export --longpi <ws>`：把本期读数、洞见读数、问题看板、方案和复测计划写成 longpi 可导入的 JSON（字段按 longpi 现有
   `save_self_measurement` / `draft_intervention_plan` / `set_followup` 的入参）。
3. 在 SKILL.md 加一段「在 dsh/longpi 中运行」：工作区位置、方法库路径默认 `~/longpi/longevity-skills`、子代理用 dsh 的 subagent。

### 3.2 longpi 侧（dsh-plugin-longpi 仓库）
1. 新工具 `run_deep_analysis`：在当前会话派发 longevity-analyst（dsh 已能发现 `~/.dsh/skills/longevity-analyst`），
   传入 Mirobody MCP 地址和会员档案；进度写入会话卡片。
2. 新工具 `import_analysis`：读取 analyst 的导出 JSON，写入 longpi 的本地 store（读数、方案草稿、随访），方案仍需会员读回确认。
3. 「健康」页新增「深度分析」tab：嵌入 report.html，单独展示问题看板和器官体检表。
4. 安装器：修 LOINC 路径；默认关闭 session-log-deepseek；增加 `--mirobody-native`；可选把 longevity-analyst 一并装到 `~/.dsh/skills`。

### 3.3 分支与依赖
- analyst：`feat/longpi-bridge`，依赖 v0.6.0（洞见层）先合入 main。
- longpi：`feat/analyst-bridge`，依赖 analyst 导出格式冻结（`la-export/1`）。
- 顺序：analyst 导出格式 → longpi import → longpi run_deep_analysis → UI tab → 端到端（李明华 + 第二个虚拟会员）。

## 4. 验收
- 虚拟会员从 longpi 页面上传 → Mirobody → 一键深度分析 → 报告和问题看板出现在 longpi → 方案读回确认 → 随访提醒生成。
- 两个不同会员各跑一遍；会话日志不外传（抓包或 dsh 日志确认）。
