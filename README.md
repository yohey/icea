
# Interface-extended CEA

インタフェース機能強化版 CEA。
頭の i は interface のつもり。

オリジナル版 CEA と同様にコマンドラインから実行できるほか，
Fortran, C++, Python からライブラリとして呼び出すことができます。

## 依存関係

あらかじめ下記がインストールされている必要があります。
* CMake 3.12 以上
* C++ コンパイラ (gcc, clang 等)
* Fortran コンパイラ (gfortran 等)

## インストール方法

```
cmake . -B _build --install-prefix=/opt/local/cea/icea-0.1.0 -DCMAKE_BUILD_TYPE=Release
cmake --build _build
cmake --install _build
```

## 動作テスト

下記により正常にビルドできたかどうかの動作テストを実行できます。
```
cmake --build _build --target test
```

ただし，コンパイラの最適化レベルが `-O2` や `-O3` の場合は計算結果が微妙に異なるためテストは失敗と判定されます。
テストを成功させたい場合は `-O0` または `-O1` を指定する必要があります。
```
cmake . -B _build --install-prefix=/opt/local/cea/icea-0.1.0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS_RELEASE="-O1"
```

デバッグしたい場合の例：
```
cmake . -B _build --install-prefix=/opt/local/cea/icea-0.1.0 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS_RELEASE="-O0 -Wall -g -fbacktrace"
```

## 使用方法

### コマンドラインから実行する方法

#### 通常モード

このモードは未実装です。

#### 再現モード

なるべくオリジナル版 CEA に近い動作をするモードです。
実行ファイルの後に `--legacy-mode` を指定し，インプットファイルは指定せずに実行します。
インプットファイル名はオリジナル版 CEA と同様に対話的な入力が求められます。

```
/opt/local/cea/icea-0.1.0/bin/cea --legacy-mode
```

### ライブラリとして呼び出す方法

#### Fortran

[examples/fortran](examples/fortran) 参照。

#### C++

[examples/c++](examples/c++) 参照。

#### Python

[examples/python](examples/python) 参照。


## 免責事項

* 本ソフトウェアは，本家 NASA CEA と同一の計算結果を保証するものではありません。
* 本ソフトウェアは使用者の自己責任において使用されるものとします。本ソフトウェアの使用に関連して発生した損害等については一切の責任を負いません。


## ライセンス

本ソフトウェアのライセンスについては検討中です (すでに NASA CEA を改造したソフトウェアは数多く公開されていますが，ものによって GPL だったり MIT ライセンスだったり様々であり，どうすればいいかわかりません)。
