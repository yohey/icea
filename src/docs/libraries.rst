物性値ライブラリについて
========================

CEA の計算実行には，thermo.lib および trans.lib という二つの物性値ライブラリを読み込む必要があります。これらについては下記の優先順位に従ってファイルを探索し，最初に発見されたファイルを読み込みます。


【優先順位 1】明示的に指定
--------------------------

物性値ライブラリのファイルパスをユーザが明示的に指定する方法です。

コマンド実行時に読み込むライブラリを指定する例::

  cea --thermo-lib <filepath> --trans-lib <filepath>


C++ コード中で読み込むライブラリを指定する例::

  prob.set_problem(mode = "rocket", name = "New Problem", thermo_lib = "<filepath>", trans_lib = "<filepath>")

明示的に指定したにもかかわらず，そのファイルが見つからない・読み込み権限がない等であった場合はエラーとなり，優先順位 2〜4 の自動探索は実施しません。


【優先順位 2】カレントディレクトリ
----------------------------------

本家 CEA と同じ挙動です。実行時のカレントディレクトリに thermo.lib と trans.lib があれば，それを読み込みます。


【優先順位 3】環境変数
----------------------

環境変数 ``CEA_ROOT`` が設定されていた場合，その階層下の ``${CEA_ROOT}/lib/cea/thermo.lib`` および ``${CEA_ROOT}/lib/cea/trans.lib`` を読み込みます。

たとえば ``CEA_ROOT`` が ``/home/user/work/cea`` に設定されていた場合は，下記のファイルが候補になります。

* ``/home/user/work/cea/lib/cea/thermo.lib``
* ``/home/user/work/cea/lib/cea/trans.lib``


【優先順位 4】インストール場所
------------------------------

当該 CEA をインストールした場所，すなわちビルド時に cmake に対して指定した場所の階層下にあるファイルを読み込みます。

たとえば，ビルド時に ``cmake . -B _build --install-prefix=/opt/local/cea/icea-0.1.0`` としていた場合は，下記のファイルが候補になります。

* ``/opt/local/cea/icea-0.1.0/lib/cea/thermo.lib``
* ``/opt/local/cea/icea-0.1.0/lib/cea/trans.lib``
