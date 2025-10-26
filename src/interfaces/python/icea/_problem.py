# -*- coding: utf-8 -*-

from ctypes import POINTER, Structure, byref, cast, c_bool, c_char_p, c_double, c_size_t, c_void_p
from typing import Iterable

from ._findlib import _libcea


class FFI_C_Ptr_Array(Structure):
    _fields_ = [('addr', c_void_p),
                ('size', c_size_t)]


class CEA_Problem:

    def __init__(self):
        _libcea.ffi_cea_new_problem.restype = c_void_p
        self._ffi = c_void_p(_libcea.ffi_cea_new_problem())

    def __del__(self):
        _libcea.ffi_cea_del_problem.argtypes = [c_void_p]
        _libcea.ffi_cea_del_problem(self._ffi)

    def set_problem(self, mode, name = None, mole_ratios = None, equilibrium = None, ions = None,
                    frozen = None, frozen_at_throat = None, incident = None, reflected = None, thermo_lib = None, trans_lib = None):
        """問題の基本設定を行います。

        Args:
            mode (str): 計算モードを設定します。有効な選択肢は ``tp``, ``pt``, ``hp``, ``ph``, ``uv``, ``vu``, ``tv``, ``vt``, ``deton``, ``shock``, ``rocket`` です。
            name (str): ケース名称を半角 15 文字以内で設定します。
            mole_ratios (bool): 混合比をモル比で定義する場合は ``True`` にします。デフォルトでは質量比で定義されます。
            equilibrium (bool): 平衡状態を仮定する場合は ``True`` にします。
            ions (bool): 燃焼過程でイオン化を考慮する場合は ``True`` にします。
            frozen (bool): ``rocket`` モード専用のオプションです。凍結流を仮定する場合は ``True`` にします。
            frozen_at_throat (bool): ``rocket`` モード専用，かつ ``frozen`` が ``True`` の場合のみ有効です。凍結位置をスロートにする場合は ``True`` にします。デフォルトの凍結位置はチャンバです。
            incident (bool): ``shock`` モード専用のオプションです。入射衝撃波を計算する場合は ``True`` にします。
            reflected (bool): ``shock`` モード専用のオプションです。反射衝撃波を計算する場合は ``True`` にします。
            thermo_lib (str): 物性値ファイル thermo.lib のを明示的に指定したい場合にフルパスを設定します。デフォルトでは自動探索します。
            trans_lib (str): 物性値ファイル trans.lib を明示的に指定したい場合にフルパスを設定します。デフォルトでは自動探索します。

        Returns:
            なし
        """
        _libcea.ffi_cea_set_problem.argtypes = [c_void_p, c_char_p, c_char_p, POINTER(c_bool), POINTER(c_bool), POINTER(c_bool),
                                                POINTER(c_bool), POINTER(c_bool), POINTER(c_bool), POINTER(c_bool), c_char_p, c_char_p]
        _libcea.ffi_cea_set_problem(self._ffi, mode.encode(), _c_char_or_None(name), _c_bool_or_None(mole_ratios), _c_bool_or_None(equilibrium), _c_bool_or_None(ions),
                                    _c_bool_or_None(frozen), _c_bool_or_None(frozen_at_throat), _c_bool_or_None(incident), _c_bool_or_None(reflected),
                                    _c_char_or_None(thermo_lib), _c_char_or_None(trans_lib))

    def set_output_options(self, SI = True, debug_points = None, mass_fractions = False, short = False, trace_tol = None, transport = False, plot = None):
        """結果出力に関する設定を行います。

        Args:
            SI (bool): 結果を SI（国際単位系）で出力します。
            debug_points (Iterable[int]): デバッグ用に冗長な出力をさせる点のインデックスを指定します。
            mass_fractions (bool): 生成物の組成情報を質量比で出力します。デフォルトではモル比で出力されます。
            short (bool): 計算の途中経過などの出力を抑制し，出力される情報量を減らします。
            trace_tol (float): モル比または質量比がこの値以上の物質のみが計算結果として出力されます。デフォルト挙動は ``shock`` モードでは 5 × 10\\ :sup:`-9`，それ以外のモードでは 5 × 10\\ :sup:`-6` となります。
            transport (bool): 出力情報に輸送関連パラメータを含めます。
            plot (Iterable[str]): plt ファイルに出力するパラメータ名を指定します。有効なキーワードは `NASA-RP-1331 <https://ntrs.nasa.gov/citations/19960044559>`_ の 2.5.4 項を参照してください。

        Returns:
            なし

        Example:
            >>> prob.set_output_options(SI = True, transport = True)
            >>> prob.set_output_options(SI = True, mass_fractions = True, plot = ['aeat', 't', 'p', 'ivac', 'isp', 'mach', 'cf'])
            >>> prob.set_output_options(SI = False, trace_tol = 1e-10)
            >>> prob.set_output_options(SI = True, debug_points = [4])
        """
        _libcea.ffi_cea_set_output_options.argtypes = [c_void_p, POINTER(c_bool), POINTER(FFI_C_Ptr_Array), POINTER(c_bool), POINTER(c_bool),
                                                       POINTER(c_double), POINTER(c_bool), c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_output_options(self._ffi, _c_bool_or_None(SI), _c_size_t_array_or_None(debug_points), _c_bool_or_None(mass_fractions),
                                           _c_bool_or_None(short), _c_double_or_None(trace_tol), _c_bool_or_None(transport), *_c_char_array_or_None(plot))

    def set_chamber_pressures(self, pressures, unit = 'MPa'):
        """チャンバ圧力を設定します。

        Args:
            pressures (Iterable[float]): 圧力の値を設定します。
            unit (str): 圧力の単位を指定します。有効な選択肢は ``MPa``，``kPa``，``Pa``，``atm``，``bar``，``psi``，``legacy-psi`` です。``psi`` は 6894.76 Pa（現代の定義）であり，``legacy-psi`` は 6894.73 Pa（オリジナル CEA に合わせたもの）です。

        Returns:
            なし

        Example:
            >>> prob.set_chamber_pressures([1.0, 2.0, 5.0], 'MPa')
        """
        _libcea.ffi_cea_set_chamber_pressures.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array), c_char_p]
        _libcea.ffi_cea_set_chamber_pressures(self._ffi, _c_double_array(pressures), _c_char_or_None(unit))

    def set_chamber_temperatures(self, temperatures, unit = 'K'):
        """チャンバ温度を設定します。

        Args:
            temperatures (Iterable[float]): 温度の値を設定します。
            unit (str): 温度の単位を設定します。現時点では有効な選択肢は ``K`` のみです。

        Returns:
            なし

        Example:
            >>> prob.set_chamber_temperatures([300.0, 500.0, 1000.0])
        """
        _libcea.ffi_cea_set_chamber_temperatures.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array), c_char_p]
        _libcea.ffi_cea_set_chamber_temperatures(self._ffi, _c_double_array(temperatures), _c_char_or_None(unit))

    def set_chamber_densities(self, densities, unit = 'kg/m^3'):
        """``tv``, ``uv``, ``sv`` モードにおいて，比体積のかわりに密度を指定する場合に使用します。

        Args:
            densities (Iterable[float]): 密度の値を設定します。
            unit (str): 密度の単位を設定します。有効な選択肢は ``kg/m^3`` または ``g/cc`` です。

        Returns:
            なし

        Example:
            >>> prob.set_chamber_densities([9.1864e-5, 8.0877e-6, 6.6054e-7], unit = 'g/cc')
        """
        _libcea.ffi_cea_set_chamber_densities.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array), c_char_p]
        _libcea.ffi_cea_set_chamber_densities(self._ffi, _c_double_array(densities), _c_char_or_None(unit))

    def set_chamber_internal_energy(self, uByR):
        """``uv`` モードにおいて，内部エネルギを指定する場合に使用します。
        内部エネルギの値は気体定数 R0 = 8.31451 [J/(K·mol)] で無次元化した値として入力します
        （現代の定義 R0 = 8.31446261815324 [J/(K·mol)] とは異なることに注意）。

        Args:
            uByR (float): 内部エネルギを気体定数で無次元化した値を設定します。

        Returns:
            なし

        Example:
            >>> prob.set_chamber_internal_energy(-45.1343)
        """
        _libcea.ffi_cea_set_chamber_internal_energy.argtypes = [c_void_p, POINTER(c_double)]
        _libcea.ffi_cea_set_chamber_internal_energy(self._ffi, _c_double_or_None(uByR))

    def set_mixture_ratios(self, ratios, type_ = None):
        """混合比を設定します。

        Args:
            ratios (Iterable[float]): 混合比の値を設定します。
            type_ (str): 混合比のタイプを指定します。デフォルトでは O/F 質量比となります。有効な選択肢は ``eq.ratio`` または ``%fuel`` です。

        Returns:
            なし

        Example:
            >>> prob.set_mixture_ratios([1.0, 1.5], type_ = 'eq.ratio')
        """
        _libcea.ffi_cea_set_mixture_ratios.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array), c_char_p]
        _libcea.ffi_cea_set_mixture_ratios(self._ffi, _c_double_array(ratios), _c_char_or_None(type_))

    def set_pressure_ratios(self, ratios):
        """圧力比を設定します。

        Args:
            ratios (Iterable[float]): 圧力比の値を設定します。

        Returns:
            なし

        Example:
            >>> prob.set_pressure_ratios([10.0, 100.0, 1000.0])
        """
        _libcea.ffi_cea_set_pressure_ratios.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_pressure_ratios(self._ffi, _c_double_array(ratios))

    def set_subsonic_area_ratios(self, ratios):
        """亜音速の面積比を設定します。

        Args:
            ratios (Iterable[float]): 面積比の値を設定します。

        Returns:
            なし

        Example:
            >>> prob.set_subsonic_area_ratios([1.58])
        """
        _libcea.ffi_cea_set_subsonic_area_ratios.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_subsonic_area_ratios(self._ffi, _c_double_array(ratios))

    def set_supersonic_area_ratios(self, ratios):
        """超音速の面積比を設定します。

        Args:
            ratios (Iterable[float]): 面積比の値を設定します。

        Returns:
            なし

        Example:
            >>> prob.set_supersonic_area_ratios([25.0, 50.0, 75.0])
        """
        _libcea.ffi_cea_set_supersonic_area_ratios.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_supersonic_area_ratios(self._ffi, _c_double_array(ratios))

    def set_finite_area_combustor(self, contraction_ratio = None, mass_flow_ratio = None):
        """有限面積の燃焼室を仮定する場合に指定します。
        引数として ``contraction_ratio`` または ``mass_flow_ratio`` どちらか一方を必ず指定する必要があります。

        Args:
            contraction_ratio (float): 収縮比（燃焼室面積／スロート面積）を指定します。
            mass_flow_ratio (float): 面積あたりの質量流量比 [(kg/s)/m^2] を指定します。

        Returns:
            なし

        Example:
            >>> prob.set_finite_area_combustor(contraction_ratio = 1.58)
        """
        _libcea.ffi_cea_set_finite_area_combustor.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
        _libcea.ffi_cea_set_finite_area_combustor(self._ffi, _c_double_or_None(contraction_ratio), _c_double_or_None(mass_flow_ratio))

    def set_initial_velocities(self, velocities, is_mach = False):
        """``shock`` モードにおいて，初期速度を指定する際に使用します。

        Args:
            velocities (Iterable[float]): 速度の値を設定します。
            is_mach (bool): ``True`` にすると指定した値をマッハ数とみなします。

        Returns:
            なし

        Example:
            >>> prob.set_initial_velocities([1000, 1100, 1200, 1250, 1300, 1350, 1400])
        """
        _libcea.ffi_cea_set_initial_velocities.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array), POINTER(c_bool)]
        _libcea.ffi_cea_set_initial_velocities(self._ffi, _c_double_array(velocities), _c_bool_or_None(is_mach))

    def add_reactant(self, type_, name, formula = None, ratio = None, T = None, rho = None, h = None, u = None,
                     T_unit = 'K', rho_unit = 'kg/m^3', h_unit = 'J/mol', u_unit = 'J/mol'):
        """新しい反応物質を追加します。``h`` と ``u`` はどちらか片方のみ指定できます。

        Args:
            type_ (str): 物質の種類を指定します。有効な選択肢は ``fuel``, ``oxyd``, ``name`` です。
            name (str): 物質名称を指定します。
            formula (str): 物質の組成をスペース区切りで指定します。たとえば CH\\ :sub:`4` ならば ``"C 1.0 H 4.0"`` とします。基本的に ``name`` に指定した名称が thermo.lib 内の名称と一致している場合は省略可能ですが，``h``, ``u`` を指定する場合は ``formula`` が必須です（今後改善予定）。
            ratio (float): 物質の存在比率を指定します。``set_problem`` で ``mole_ratios`` を ``True`` にした場合はモル比，そうでない場合は質量比として解釈されます。値は自動的に正規化されます (たとえば 0.4, 0.6 のように指定しても 40, 60 のように指定しても同じ結果となります)。
            T (float): 物質の温度を指定します。
            rho (float): 物質の密度を指定します。
            h (float): 物質のエンタルピを指定します。
            u (float): 物質の内部エネルギを指定します。
            T_unit (str): 温度の単位を指定します。有効な選択肢は現時点では ``K`` のみです。
            rho_unit (str): 密度の単位を指定します。有効な選択肢は ``kg/m^3``, ``g/cm^3``, ``g/cc`` です。
            h_unit (str): エンタルピの単位を指定します。有効な選択肢は ``J/mol``, ``kJ/mol``, ``cal/mol`` です。
            u_unit (str): 内部エネルギの単位を指定します。有効な選択肢は ``J/mol``, ``kJ/mol``, ``cal/mol`` です。

        Returns:
            なし

        Example:
            >>> prob.add_reactant('fuel', 'H2(L)', ratio = 100, T = 20.27)
            >>> prob.add_reactant('oxyd', 'O2(L)', ratio = 100, T = 90.17)
        """
        _libcea.ffi_cea_add_reactant.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p, POINTER(c_double),
                                                 POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                                 c_char_p, c_char_p, c_char_p, c_char_p]
        _libcea.ffi_cea_add_reactant(self._ffi, type_.encode(), name.encode(), _c_char_or_None(formula), _c_double_or_None(ratio),
                                     _c_double_or_None(T), _c_double_or_None(rho), _c_double_or_None(h), _c_double_or_None(u),
                                     _c_char_or_None(T_unit), _c_char_or_None(rho_unit), _c_char_or_None(h_unit), _c_char_or_None(u_unit))

    def set_reactant(self, index, ratio = None, T = None, rho = None, h = None, u = None,
                     T_unit = 'K', rho_unit = 'kg/m^3', h_unit = 'J/mol', u_unit = 'J/mol'):
        """既存の反応物質の設定を変更します。

        Args:
            index (int): 変更したい物質のインデックス（0 始まりで ``add_reactant`` した順）を指定します。
            ratio (float): 物質の存在比率を指定します。
            T (float): 物質の温度を指定します。
            rho (float): 物質の密度を指定します。
            h (float): 物質のエンタルピを指定します。
            u (float): 物質の内部エネルギを指定します。
            T_unit (str): 温度の単位を指定します。
            rho_unit (str): 密度の単位を指定します。
            h_unit (str): エンタルピの単位を指定します。
            u_unit (str): 内部エネルギの単位を指定します。

        Returns:
            なし
        """
        _libcea.ffi_cea_set_reactant.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_double),
                                                 POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                                 c_char_p, c_char_p, c_char_p, c_char_p]
        _libcea.ffi_cea_set_reactant(self._ffi, byref(c_size_t(index + 1)), _c_double_or_None(ratio),
                                     _c_double_or_None(T), _c_double_or_None(rho), _c_double_or_None(h), _c_double_or_None(u),
                                     _c_char_or_None(T_unit), _c_char_or_None(rho_unit), _c_char_or_None(h_unit), _c_char_or_None(u_unit))

    def set_insert_species(self, species):
        """最初の計算点から考慮に追加する物質を指定します。

        Args:
            species (Iterable[str]): 考慮に追加する物質を指定します。thermo.lib 内の名称と完全に一致させる必要があります。

        Returns:
            なし

        Example:
            >>> prob.set_insert_species(['BeO(L)'])
        """
        _libcea.ffi_cea_set_insert_species.argtypes = [c_void_p, c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_insert_species(self._ffi, *_c_char_array_or_None(species))

    def set_omit_species(self, species):
        """燃焼過程において無視する物質を設定します。

        Args:
            species (Iterable[str]): 無視する物質を指定します。thermo.lib 内の名称と完全に一致させる必要があります。

        Returns:
            なし

        Example:
            >>> prob.set_omit_species(['C(gr)'])
        """
        _libcea.ffi_cea_set_omit_species.argtypes = [c_void_p, c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_omit_species(self._ffi, *_c_char_array_or_None(species))

    def set_only_species(self, species):
        """計算に考慮する物質を限定する場合に使用します。

        Args:
            species (Iterable[str]): 計算に考慮する物質を指定します。thermo.lib 内の名称と完全に一致させる必要があります。

        Returns:
            なし

        Example:
            >>> prob.set_only_species(['Ar', 'CO', 'CO2', 'H2', 'H2O', 'HNO', 'HO2', 'NH', 'NO', 'N2', 'O2', 'OH'])
        """
        _libcea.ffi_cea_set_only_species.argtypes = [c_void_p, c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_only_species(self._ffi, *_c_char_array_or_None(species))

    def set_legacy_mode(self, legacy_mode = None):
        """なるべく本家 CEA に近い動作を再現したい場合に指定します。本家 CEA と同じ結果を保証するわけではありません。

        Args:
            legacy_mode (bool): ``True`` にすると再現モードとなります。

        Returns:
            なし

        Example:
            >>> prob.set_legacy_mode(True)
        """
        _libcea.ffi_cea_set_legacy_mode.argtypes = [c_void_p, POINTER(c_bool)]
        _libcea.ffi_cea_set_legacy_mode(self._ffi, _c_bool_or_None(legacy_mode))

    def run(self, out_filename = None, plt_filename = None):
        """計算を実行します。

        Args:
            out_filename (str): アウトプットファイル名を指定します。省略した場合は，
                                ``legacy_mode = True`` ならば標準出力に結果を出力し，
                                ``legacy_mode = False`` ならば何も出力せずに計算だけが実行されます。
            plt_filename (str): プロットファイル名を指定します。省略した場合は何も出力しません。

        Returns:
            なし

        Example:
            >>> prob.run(out_filename = 'example-12.out', plt_filename = 'example-12.plt')
        """
        _libcea.ffi_cea_run.argtypes = [c_void_p, c_char_p, c_char_p]
        _libcea.ffi_cea_run(self._ffi, _c_char_or_None(out_filename), _c_char_or_None(plt_filename))

    def get_pressure(self, iOF, ipt):
        """計算結果の圧力の値を取得します。``.run()`` の実行より後で使用します。

        Args:
            iOF (int): 0 から始まる混合比インデックスを指定します。
            ipt (int): 0 から始まる計算点インデックスを指定します。

        Returns:
            float: 圧力 [MPa]

        Example:
            >>> P = prob.get_pressure(0, 2)
        """
        _libcea.ffi_cea_get_pressure.restype = c_double
        _libcea.ffi_cea_get_pressure.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_pressure(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_temperature(self, iOF, ipt):
        """計算結果の温度の値を取得します。``.run()`` の実行より後で使用します。

        Args:
            iOF (int): 0 から始まる混合比インデックスを指定します。
            ipt (int): 0 から始まる計算点インデックスを指定します。

        Returns:
            float: 温度 [K]

        Example:
            >>> T = prob.get_temperature(0, 2)
        """
        _libcea.ffi_cea_get_temperature.restype = c_double
        _libcea.ffi_cea_get_temperature.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_temperature(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_chamber_temperature(self):
        """``rocket`` モードにおいて，計算結果のチャンバ温度の値を取得します。``.run()`` の実行より後で使用します。
        ``FAC`` (Finite Area Combustor) モードが有効の場合は ``get_temperature(0, 1)`` と同じ，
        それ以外の場合は ``get_temperature(0, 0)`` と同じです。

        Args:
            なし

        Returns:
            float: 温度 [K]

        Example:
            >>> Tc = prob.get_chamber_temperature()
        """
        _libcea.ffi_cea_get_chamber_temperature.restype = c_double
        _libcea.ffi_cea_get_chamber_temperature.argtypes = [c_void_p]
        return _libcea.ffi_cea_get_chamber_temperature(self._ffi)

    def get_molecular_weight(self, iOF, ipt):
        """計算結果のモル質量の値を取得します。``.run()`` の実行より後で使用します。

        Args:
            iOF (int): 0 から始まる混合比インデックスを指定します。
            ipt (int): 0 から始まる計算点インデックスを指定します。

        Returns:
            float: モル質量 [g/mol]

        Example:
            >>> M = prob.get_molecular_weight(0, 2)
        """
        _libcea.ffi_cea_get_molecular_weight.restype = c_double
        _libcea.ffi_cea_get_molecular_weight.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_molecular_weight(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_specific_heat(self, iOF, ipt):
        """計算結果の比熱の値を取得します。``.run()`` の実行より後で使用します。

        Args:
            iOF (int): 0 から始まる混合比インデックスを指定します。
            ipt (int): 0 から始まる計算点インデックスを指定します。

        Returns:
            float: 比熱 [kJ/(kg·K)] または [cal/(g·K)] （単位は ``set_output_options`` の ``SI`` により変化）

        Example:
            >>> cp = prob.get_specific_heat(0, 2)
        """
        _libcea.ffi_cea_get_specific_heat.restype = c_double
        _libcea.ffi_cea_get_specific_heat.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_specific_heat(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_specific_heat_ratio(self, iOF, ipt):
        """計算結果の比熱比の値を取得します。``.run()`` の実行より後で使用します。

        Args:
            iOF (int): 0 から始まる混合比インデックスを指定します。
            ipt (int): 0 から始まる計算点インデックスを指定します。

        Returns:
            float: 比熱比 [-]

        Example:
            >>> gamma = prob.get_specific_heat_ratio(0, 2)
        """
        _libcea.ffi_cea_get_specific_heat_ratio.restype = c_double
        _libcea.ffi_cea_get_specific_heat_ratio.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_specific_heat_ratio(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_characteristic_velocity(self):
        """``rocket`` モードにおいて，計算結果の特性排気速度の値を取得します。``.run()`` の実行より後で使用します。

        Args:
            なし

        Returns:
            float: 特性排気速度 [m/s]

        Example:
            >>> cstar = prob.get_characteristic_velocity()
        """
        _libcea.ffi_cea_get_characteristic_velocity.restype = c_double
        _libcea.ffi_cea_get_characteristic_velocity.argtypes = [c_void_p]
        return _libcea.ffi_cea_get_characteristic_velocity(self._ffi)

    def get_specific_impulse(self, iOF, ipt, vacuum = False):
        """``rocket`` モードにおいて，計算結果の比推力の値を取得します。``.run()`` の実行より後で使用します。

        Args:
            iOF (int): 0 から始まる混合比インデックスを指定します。
            ipt (int): 0 から始まる計算点インデックスを指定します。
            vacuum (bool): ``True`` の場合は真空中比推力を返します。

        Returns:
            float: 比推力 [s]

        Example:
            >>> Isp = prob.get_specific_impulse(0, 2)
        """
        _libcea.ffi_cea_get_specific_impulse.restype = c_double
        _libcea.ffi_cea_get_specific_impulse.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t), POINTER(c_bool)]
        return _libcea.ffi_cea_get_specific_impulse(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)), _c_bool_or_None(vacuum))

    def calc_frozen_exhaust(self, T):
        """凍結流として指定の温度になった場合の物性値を取得します。``.run()`` の実行より後で使用します。

        Args:
            T (float): 温度 [K]

        Returns:
            tuple[float, float, float, float, float]:
            比熱 [kJ/(kg·K)]，比熱比 [-]，粘性係数 [µPa·s]，熱伝導率 [W/(m·K)]，プラントル数 [-]

        Example:
            >>> cp, gamma, mu, k, Pr = prob.calc_frozen_exhaust(650.)
        """
        _libcea.ffi_cea_calc_frozen_exhaust.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        _cp, _gamma, _mu, _k, _Pr = c_double(), c_double(), c_double(), c_double(), c_double()
        _libcea.ffi_cea_calc_frozen_exhaust(self._ffi, _c_double_or_None(T), byref(_cp), byref(_gamma), byref(_mu), byref(_k), byref(_Pr))
        return _cp.value, _gamma.value, _mu.value, _k.value, _Pr.value

    def get_thermo_reference_properties(self, name: str):
        """指定した物質の参照状態量を取得します。

        Args:
            name (str): 物質名

        Returns:
            tuple[float, float, float]: モル質量 [g/mol]，参照温度 [K]，参照エンタルピ [kJ/kg]

        Example:
            >>> M, T_ref, h_ref = prob.get_thermo_reference_properties('O2(L)')
        """
        _libcea.ffi_cea_get_thermo_reference_properties.argtypes = [c_void_p, c_char_p, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        _M, _T_ref, _h0_ref = c_double(), c_double(), c_double()
        _libcea.ffi_cea_get_thermo_reference_properties(self._ffi, name.encode(), byref(_M), byref(_T_ref), byref(_h0_ref))
        return _M.value, _T_ref.value, _h0_ref.value


    def write_debug_output(self, filename: str):
        """デバッグ用の情報をファイルに出力します。

        Args:
            filename (str): 出力ファイル名

        Returns:
            なし
        """
        _libcea.ffi_cea_write_debug_output.argtypes = [c_void_p, c_char_p]
        _libcea.ffi_cea_write_debug_output(self._ffi, filename.encode())


def _c_bool_or_None(b: bool | None):
    if b is None:
        return None
    else:
        return byref(c_bool(b))

def _c_double_or_None(x: float | None):
    if x is None:
        return None
    else:
        return byref(c_double(x))

def _c_char_or_None(s: str | None):
    if s is None:
        return None
    else:
        return s.encode()

def _c_double_array(a: Iterable[float]):
    _a = (c_double * len(a))(*a)
    ptr = FFI_C_Ptr_Array(addr = cast(_a, c_void_p), size = len(a))
    return byref(ptr)

def _c_size_t_array_or_None(a: Iterable[int] | None):
    if a is None:
        return None
    else:
        _a = (c_size_t * len(a))(*(i + 1 for i in a))
        ptr = FFI_C_Ptr_Array(addr = cast(_a, c_void_p), size = len(a))
        return byref(ptr)

def _c_char_array_or_None(a: Iterable[str] | None):
    if a is None:
        return None, None
    else:
        char_array_ptr = (c_char_p * len(a))(*[s.encode() for s in a])
        char_length_array = (c_size_t * len(a))(*[len(s) for s in a])
        char_length_array_ptr = FFI_C_Ptr_Array(addr = cast(char_length_array, c_void_p), size = len(a))
        return byref(char_array_ptr), byref(char_length_array_ptr)
