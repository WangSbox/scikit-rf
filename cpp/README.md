# skrf C++14 minimal port

这是一个最小化的 C++14 项目骨架，用于把 `scikit-rf` 的部分功能迁移到 C++。

包含：
- 一个简单的 `Network` 结构
- 基本的 Touchstone (`.sNp`) 文件解析器（有限的、用于示例）
- `CMakeLists.txt` 和 `src/main.cpp` 示例

构建与运行：

```
mkdir build && cd build
cmake ..
cmake --build .
./skrf_sample ../skrf/data/ntwk1.s2p
```

说明：当前解析器非常简化，只支持每行完整的数字记录（频率后跟 2* N*N 个实数列），用于快速原型和后续扩展。
