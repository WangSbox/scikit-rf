# skrf 从 Python 转到 C++14 — README

## 项目简介

将 Python 版的 **scikit-rf** 迁移为 **C++14** 的工程化实现，目标是把常用功能（Touchstone I/O、S/Z/Y/ABCD 参数变换、网络级联、插值、基本校准工具等）全部用 C++ 实现，便于嵌入到高性能测量与仿真流程中。

该 README 说明了工程目标、必要依赖、代码组织、API 对应表、构建/测试/文档流程、性能/精度注意事项及贡献规范等。源码目录为sfrf文件夹。

## 目标

* 功能等价：把 Python `scikit-rf` 的常用功能（在你的代码库中实际使用到的模块）均转换为 C++14 实现。所有原来用 Python 实现的逻辑应当以相应的 C++ 类/函数替代。第三方库允许使用，但必须是现实存在并可获取的开源或商业库。
* 可用性：提供清晰的 CMake 构建、单元测试和 API 文档。
* 可嵌入：API 设计适于在测量控制/仿真管线中被直接调用（无 Python 依赖）。
* 移除绘图等组件特性。

## 强制性依赖（建议版本范围）

* 核心线性代数：entity ["organization","Eigen","linear algebra library"] (建议 >= 3.3，header-only)
* 常用 C++ 工具/组件：entity ["organization","Boost","c++ libraries"]（用于可选的字符串、filesystem、optional 等；建议 >= 1.65）
* 构建系统：entity ["organization","CMake","build system"] (建议 >= 3.10)
* 单元测试（任选）：建议使用 entity ["organization","Catch2","cpp testing framework"] 或 entity ["organization","GoogleTest","cpp testing framework"]
* 文档生成（可选但强烈建议）：entity ["organization","Doxygen","documentation generator"]
* 包管理（可选）：推荐 entity ["organization","Conan","cpp package manager"] 或 entity ["organization","vcpkg","microsoft package manager"] 以便 CI 中统一依赖管理。
* 源代码托管/Issue：项目建议托管于 entity ["organization","GitHub","code hosting platform"]。
* （背景）原始参考实现：entity ["organization","scikit-rf","python package"] — 用于功能/行为对照与测试矢量生成。

> 说明：上面列出的库都是可获取的真实项目；在 `CMake` 中以 `find_package` / 包管理器方式引用。

## 代码组织-源码
```
skrf/
```

## 代码组织cpp（建议目录结构）

```
cpp/
├─ CMakeLists.txt
├─ src/
│  ├─ skrf/                # 公共 API (namespace skrf)
│  │  ├─ frequency.hpp
│  │  ├─ network.hpp
│  │  ├─ touchstone.hpp
│  │  ├─ sparams.hpp
│  │  └─ io/               # touchstone 解析器与写入器
│  │     └─ touchstone.cpp
|  |  |_ ...
|  | 
│  └─ tools/               # 变换 / utils
│     ├─ transforms.cpp
│     └─ interpolation.cpp
├─ include/                # （若需要）对外头文件分发
├─ tests/
│  ├─ unit/                # 单元测试
│  └─ data/                # 用于测试的 touchstone 文件（从 scikit-rf 导出）
├─ docs/                   # Doxygen 配置与扩展文档
└─ examples/
   └─ read_touchstone_example.cpp
```

## API 设计要点（从 Python → C++ 的映射建议）

* `Frequency`（Python：`skrf.frequency.Frequency`）

  * C++：`skrf::Frequency` 类，持有 `std::vector<double>`（Hz）。提供插值辅助与边界/索引操作。
* `Network`（Python：`skrf.network.Network`）

  * C++：`skrf::Network`，持有 `size_t ports`、`Frequency freq` 与 `std::vector<Eigen::MatrixXcd> s_params`。
  * 方法：`s_at_index(size_t)`, `s_interp(double freq_hz)`, `to_z()`, `to_y()`, `to_abcd()` 等。
* Touchstone I/O

  * C++：`skrf::io::Touchstone` 类或 free functions：`read_touchstone(path)` 返回 `Network`；`write_touchstone(path, Network)` 写出。
  * 支持格式：RI / MA / DB 与常见频率单位（Hz/kHz/MHz/GHz）。
* 变换与网络算子

  * 提供 `s_to_z`, `z_to_s`, `s_to_abcd`, `abcd_to_s` 等函数（使用 `Eigen::Matrix2cd` / `Eigen::MatrixXcd` 作为矩阵类型）。
* 例子/兼容层

  * 提供用于与 Python 版结果对比的序列化（例如 CSV 或二进制测试向量），便于自动化回归测试。

## 精度与数值注意事项

* 默认使用 `double` 与 `std::complex<double>`，并用 `Eigen::MatrixXcd` 存储复矩阵。
* 矩阵运算（行列数可变）采用 `Eigen` 优化路径；对 2×2 常用转换可以专门实现高效版本以避免动态分配。
* 插值：建议至少支持线性插值（简单可靠），并提供可选的 cubic spline（可用 `Boost.Math` 或自实现的自然样条）。

## 构建与示例（本地）

1. 安装依赖（示例，按平台调整）：

   * 安装 entity ["organization","Eigen","linear algebra library"]（系统包或手动下载）
   * 安装 entity ["organization","CMake","build system"]
   * 可选：通过 entity ["organization","Conan","cpp package manager"] 或 entity ["organization","vcpkg","microsoft package manager"] 获取 Boost / Catch2 等

2. 构建：

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -- -j$(nproc)
```

3. 运行示例：

```bash
./examples/read_touchstone_example path/to/file.s2p
```

## 单元测试与回归测试

* 为每一个被移植的 Python 功能编写等价的 C++ 单元测试；测试数据（Touchstone 文件、预期 S 参数矩阵）最好从原始 `scikit-rf` 生成并保存在 `tests/data` 中。
* 推荐测试框架：选择 `Catch2` 或 `GoogleTest`，并在 `CMake` 中提供开关 `-DENABLE_TESTS=ON`。
* 在 CI 中运行如下测试用例：

  * I/O 正确性（读/写不改变数据格式）
  * S↔Z↔Y↔ABCD 相互转换的一致性
  * 2-port 级联运算与已知例子的对比
  * 插值结果与 Python 参考实现误差在可接受范围（例如 1e-12 绝对误差或相对误差）

## 文档

* 使用 entity ["organization","Doxygen","documentation generator"] 生成 API 文档，文档注释遵循 Doxygen 语法。
* 在 `docs/` 中放置使用示例（包括如何读取 Touchstone、获取插值 S、做网络级联与转换）的 C++ 代码片段。

## 性能与优化建议

* 对小矩阵（2×2、4×4）专门实现固定大小（`Eigen::Matrix2cd`/`Matrix4cd`）的快速路径，减少动态内存分配。
* 使用 `-O3` 与 `-march=native`（在 CI/发布构建中可选）以提高数值性能。
* 在热点处考虑 SIMD/向量化或批量处理多个频点的并行化（OpenMP 或线程池），并确保数值一致性。

## 兼容性与平台支持

* 目标支持平台：Linux（x86_64）、Windows（MSVC）与 macOS（Clang）。
* 要求 C++14 编译器支持（GCC >= 5.4 / Clang >= 3.8 / MSVC 2017 以上为常见基线）。

## 限制与差异说明（相较于 Python 版本）

* Python 版的动态特性（如 `pandas`/`numpy` 风格的广播、duck typing）在 C++ 中将被静态类型与明确 API 替代；一些便捷用法需用适配层（helper functions）重现。
* 一些高级功能（如 plotting、interactive notebooks）在 C++ 中不直接提供；建议通过生成数据文件并用 Python / MATLAB /其他工具绘图。

## 迁移流程建议（步骤化）

1. **功能梳理**：确定你当前 Python 代码中实际使用到的 `scikit-rf` 子集（列出模块/函数）。
2. **优先级划分**：优先实现 I/O、数据结构、基本变换与测试用例；把高阶功能（校准、拟合）列为 v2。
3. **接口设计**：为每个 Python API 设计 C++ 对等签名，写出 header（`.hpp`）草案并审查。
4. **实现与测试**：实现单元并用 Python 生成的测试向量验证数值一致性。
5. **性能调优**：基准测试并优化热点。
6. **文档与示例**：补充 Doxygen 文档与使用示例。

## 贡献指南

* Fork -> feature-branch -> Pull Request。
* PR 中必须包含单元测试/回归测试与文档（若修改 API）。
* 风格：遵循项目中的 clang-format 或一个明确的代码风格配置。

## 许可证

* 建议与原 Python 项目兼容的开源许可证（例如 MIT / BSD / Apache 2.0 等），并在 `LICENSE` 中明确标注。若原 `scikit-rf` 的许可需要保留声明，请在移植中保留必要的版权/许可信息。

## 常见问题（FAQ）

* **能否用 Python/C++ 混合实现？**
  可以，但本 README 的目标是“全部转换为 C++”。如果需要 Python 接口，建议额外提供一个薄的绑定层（`pybind11`），而不是把核心仍留在 Python 中。

* **必须使用哪些第三方库？**
  只要这些库在现实世界中存在并且可安装即可（比如上面列举的库）。项目应优先使用 header-only 或广泛可用的库以减少部署成本。

* **如何验证数值一致性？**
  使用 `scikit-rf` 在 Python 中生成一组标准测试向量（Touchstone 文件、S 矩阵 CSV 等），把这些向量加入 `tests/data`，在 C++ 测试中读取并断言误差界限。

-------------------------
已实现（C++）进度与单文件包含说明
-------------------------

当前仓库已完成并实现以下主要模块与特性（C++14，Eigen 作为线性代数后端）：

- 单文件入口：`cpp/include/skrf_cpp/skrf_cpp.h`，包含核心 API（`Network`、`NetworkEigen`、`Touchstone`、`media` 模型、`VectorFittingAdvanced`、`NetworkSet`、`Calibration` 等），可通过包含该头快速开始。
- Touchstone I/O：支持读取/写入 `.sNp`（支持 MA/DB/RI 格式与 Z0 处理）。
- Network / NetworkEigen：提供扁平与矩阵两种网络表示与互转，包含插值、重排、端口抽取、端口合并（merge_ports）、重采样等操作。
- 变换工具：实现 `s_to_z`, `z_to_s`, `s_to_abcd_2port`, `abcd_to_s_2port` 等常用转换，提供高效的 2×2 路径。
- 网络算子：级联（cascade）、并联（parallel）、块对角合并等。
- 媒体模型：已将多媒体模型合并到单个头文件 `cpp/include/skrf_cpp/media.h`，包括 `FreeSpace`, `CoaxialFull`, `Microstrip`, `CPW`, `RectangularWaveguide`, `CircularWaveguide` 等。
- 传输线工具：`MLine`（ABCD 矩阵生成）与设备级别组合工具。
- 校准（Calibration / CalibrationSet）：基本去嵌入（de-embedding）流程的实现与批量处理接口。
- 向量拟合（VectorFittingAdvanced）：实现迭代的 Vector Fitting 算法（可设置极点数、迭代次数、容差、正则化与权重），并暴露对单个 S 元素拟合的便捷接口。
- NetworkSet：批量网络容器，提供 `mean_s` / `std_s`、块对角合并、按索引级联/并联、以及在集合内合并端口的便捷方法。

移除与未实现的部分
-----------------

- 绘图与 Notebook 相关功能（Python 版的 matplotlib / notebook 交互）在 C++ 版本中已移除。若需要可由外部工具读取导出的 CSV/Touchstone 并绘图。
- 当前未添加 Python 绑定或示例 CI（可按需添加）。

快速示例：以单头文件开始
-----------------------

以下代码展示如何使用单文件包含来读取 Touchstone，并对 S21 进行向量拟合：

```cpp
#include "skrf_cpp.h"
using namespace skrf_cpp;

int main(){
    // 读取 touchstone 文件
    auto net = Touchstone::load("data/ntwk1.s2p");

    // 拟合 S21（端口索引 1->0，即 S(2,1) in 1-based）
    VectorFittingAdvanced vf;
    VectorFittingAdvanced::FitOptions opts;
    opts.npoles = 8;
    opts.max_iter = 50;
    vf.fit_network(net, /*port_row=*/1, /*port_col=*/0, opts);

    // 输出拟合值到文件（示例 runner 已集成）
    return 0;
}




