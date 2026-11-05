# 提交 conda-forge（staged-recipes）

本目录只保存提交 conda-forge 时用的 recipe。**上架成功后请以
`conda-forge/specproc-feedstock` 仓库里的 recipe 为准**，那里的版本由 conda-forge bot
在每次 PyPI 发版后自动更新，本文件不会跟着自动变。

> ⚠️ **格式很重要**：conda-forge 的 staged-recipes 已迁移到 CEP-13 的
> **`recipe.yaml`（v1 格式，rattler-build 构建）**，旧的 `meta.yaml`（v0）已弃用。
> 在同一个 `recipes/` 目录下混用两种格式会直接报错：
> `ValueError: Mixing meta.yaml and recipe.yaml recipes is not supported`。

## 前提

- 包已在 PyPI 上（`pip install specproc` 可用）
- sdist 干净、体积合理（当前 1.0.2 约 1.3 MB）
- 已声明 `license`、`home`、依赖，且 `pip check` 通过

## 提交步骤

1. **Fork** <https://github.com/conda-forge/staged-recipes>（网页右上角 Fork）

2. 在你 fork 的仓库里，新建分支（例如 `add-specproc`），把本目录的
   **`recipe/recipe.yaml`** 放到仓库的 **`recipes/specproc/recipe.yaml`**
   （注意：不是 `recipes/example/`，要新建 `recipes/specproc/` 目录）

3. **填好 sha256**：打开 <https://pypi.org/project/specproc/#files>，
   在 `specproc-<版本>.tar.gz` 那一行复制 SHA256，替换 recipe 里的
   `TODO_FILL_AFTER_PYPI_UPLOAD`。
   （也可以用 `sha256sum specproc-1.0.2.tar.gz` 自己算，两者必须一致）

4. 提 **Pull Request** 到 `conda-forge/staged-recipes` 的 `main` 分支：
   - 标题：`Add specproc`
   - 描述里勾选模板给出的 checklist（"I have read the guidelines" 等）

5. 等 CI：会自动在 Linux / macOS / Windows 上构建并跑 `test:` 段。
   失败时点开日志看 `conda-build` 的报错，改好后 push 到同一分支即可重跑。

6. Reviewer 通过后合并，conda-forge 会自动创建 **`specproc-feedstock`** 仓库，
   并在几分钟后出现在 `conda-forge` 频道。

## 上架之后

- 安装方式变成：`conda install -c conda-forge specproc`
- 之后每次你在 PyPI 发新版本，conda-forge bot 会自动向 feedstock 提 PR 升级，
  你只需要在 feedstock 的 PR 里等 CI 通过并 merge（通常是点一下）
- 记得把 `README.md` / `README_CN.md` 里的 conda 安装说明改回
  `conda install -c conda-forge specproc`

## 常见失败原因

| 现象 | 原因 | 处理 |
|---|---|---|
| `Mixing meta.yaml and recipe.yaml recipes is not supported` | 用了已弃用的 v0 格式 | 改名为 `recipe.yaml` 并换成 v1 语法（本目录文件即是） |
| `sha256 mismatch` | sdist 版本或哈希填错 | 用 PyPI 上对应版本的 SHA256 重新核对 |
| `python_min` 相关 lint 报错 | recipe 里 `python_min` 低于 conda-forge 当前下限 | 目前用 3.9；若被拒按提示调整 |
| `ModuleNotFoundError` | run 依赖漏了 | 补到 `requirements.run`（注意 conda 里 PyQt5 叫 `pyqt`） |
| `pip check` 失败 | 依赖版本冲突 | 看报错补版本约束 |
| 测试里 `specproc` 卡住/失败 | 无显示环境 | `main.py` 只在启动 GUI 时才 import PyQt5（第 578 行），`specproc --version` 不触发；测试仍以 `import` + `pip_check` 为主 |

## 依赖映射（PyPI → conda-forge）

| PyPI | conda-forge |
|---|---|
| `PyQt5` | `pyqt` |
| `matplotlib` | `matplotlib-base`（不拉 GUI 后端，`pyqt` 单独声明） |
| `numpy` / `scipy` / `astropy` | 同名 |
| ~~`configparser`~~ | 已删除（Python 3.2 起为标准库） |
