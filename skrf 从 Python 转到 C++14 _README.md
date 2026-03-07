<!--
 * @Author: WangSibo, wsb304617488@gmail.com
 * @Date: 2026-03-06 21:52:15
 * @LastEditors: WangSibo
 * @LastEditTime: 2026-03-07 17:52:08
 * @FilePath: \scikit-rf\skrf 从 Python 转到 C++14 _README.md
 * @Description: Copyright (c) 2026 by wsb304617488@gmail.com, All Rights Reserved.
 * it's ok, everything will be ok.
 * Copyright (c) 2026 by wsb304617488@gmail.com, All Rights Reserved. 
-->
# skrf 从 Python 转到 C++14 — README

## 项目简介

将 Python 版的 **scikit-rf** 迁移为 **C++14** 的工程化实现，目标是把常用功能（Touchstone I/O、S/Z/Y/ABCD 参数变换、网络级联、插值、基本校准工具等）全部用 C++ 实现，便于嵌入到高性能测量与仿真流程中。

该 README 说明了工程目标、必要依赖、代码组织、API 对应表、构建/测试/文档流程、性能/精度注意事项及贡献规范等，源码目录为sfrf文件夹。

## 目标

* 功能等价：把 Python `scikit-rf` 的常用功能均转换为 C++14 实现。所有原来用 Python 实现的逻辑应当以相应的 C++ 类/函数替代。第三方库允许使用，但必须是现实存在并可获取的开源或商业库。
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

## 单元测试与回归测试（暂不实现）

* 为每一个被移植的 Python 功能编写等价的 C++ 单元测试；测试数据（Touchstone 文件、预期 S 参数矩阵）最好从原始 `scikit-rf` 生成并保存在 `tests/data` 中。
* 推荐测试框架：选择 `Catch2` 或 `GoogleTest`，并在 `CMake` 中提供开关 `-DENABLE_TESTS=ON`。
* 在 CI 中运行如下测试用例：

  * I/O 正确性（读/写不改变数据格式）
  * S↔Z↔Y↔ABCD 相互转换的一致性
  * 2-port 级联运算与已知例子的对比
  * 插值结果与 Python 参考实现误差在可接受范围（例如 1e-12 绝对误差或相对误差）

## 性能与优化建议

* 对小矩阵（2×2、4×4）专门实现固定大小（`Eigen::Matrix2cd`/`Matrix4cd`）的快速路径，减少动态内存分配。
* 使用 `-O3` 与 `-march=native`（在 CI/发布构建中可选）以提高数值性能。
* 在热点处考虑 SIMD/向量化或批量处理多个频点的并行化（OpenMP 或线程池），并确保数值一致性。

## 兼容性与平台支持

* 目标支持平台：Linux（x86_64）、Windows（MSVC）与 macOS（Clang）。
* 要求 C++14 编译器支持（GCC >= 5.4 / Clang >= 3.8 / MSVC 2017 以上为常见基线）。

## 限制与差异说明（相较于 Python 版本）

* Python 版的动态特性（如 `pandas`/`numpy` 风格的广播、duck typing）在 C++ 中将被静态类型与明确 API 替代；一些便捷用法需用适配层重现。
* 一些高级功能（如 plotting、interactive notebooks）在 C++ 中不直接提供；建议通过生成数据文件并用 Python / MATLAB /其他工具绘图。

## 迁移流程建议（步骤化）

1. **功能梳理**：确定你当前 Python 代码中实际使用到的 `scikit-rf` 子集（列出模块/函数）。
2. **优先级划分**：优先实现 I/O、数据结构、基本变换与测试用例；把高阶功能（校准、拟合）列为 v2。
3. **接口设计**：为每个 Python API 设计 C++ 对等签名，写出 header（`.hpp`）草案并审查。
4. **实现与测试**：实现单元并用 Python 生成的测试向量验证数值一致性。
5. **性能调优（暂不实现）**：基准测试并优化热点。
6. **文档与示例（暂不实现）**：补充 Doxygen 文档与使用示例。

