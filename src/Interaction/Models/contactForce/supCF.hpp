#ifndef __supCF_hpp__
#define __supCF_hpp__

#include "symArrays.hpp"
#include "types.hpp"

namespace pFlow::cfModels {

/**
 * @brief SUP&JKR-F黏附力模型实现
 *
 * JKR-F模型使用表面能密度γ（J/m²）而非体积能密度Cd（J/m³）
 *
 * 原始粒子尺度的JKR-F法向力：
 * F_n^O = (4/3)√(R_eff*E*)δ^(3/2) - √(16πγE*)(R_eff*δ)^(3/4)
 *
 * SUP缩放关系：
 * - 半径：R_S = l × R_O
 * - 重叠：δ_S = l × δ_O
 * - 力：F_S = l² × F_O
 * - 表面能密度γ在两个尺度下保持不变
 *
 * @tparam limited 是否限制切向力（默认true）
 */

template <bool limited = true> class sup {
public:
  /**
   * @brief 接触历史存储
   */
  struct contactForceStorage {
    realx3 overlap_t_ = 0.0; 			// 切向重叠历史（缩放粒子坐标系）
  };

  /**
   * @brief SUP&JKR-F属性
   */
  struct supProperties {
    // 缩放因子
    real l_ = 1.0; // l = d_S/d_O（缩放粒子直径/原始粒子直径）

    // 材料参数（每种材料对）
    real Yeff_ = 1e7;   // 有效杨氏模量 [Pa]
    real Geff_ = 4e6;   // 有效剪切模量 [Pa]
    real gamma_ = 0.01; // 表面能密度 [J/m²]
    real en_ = 0.8;     // 法向恢复系数
    real et_ = 0.8;     // 切向恢复系数
    real mu_ = 0.5;     // 滑动摩擦系数

    INLINE_FUNCTION_HD
    supProperties() {}

    INLINE_FUNCTION_HD
    supProperties(real l, real Yeff, real Geff, real gamma, real en,
                      real et, real mu)
        : l_(l), Yeff_(Yeff), Geff_(Geff), gamma_(gamma), en_(en), et_(et),
          mu_(mu) {}

    INLINE_FUNCTION_HD
    supProperties(const supProperties &) = default;

    INLINE_FUNCTION_HD
    supProperties &operator=(const supProperties &) = default;

    INLINE_FUNCTION_HD
    ~supProperties() = default;
  };

protected:
  using PropertiesArrayType = symArray<supProperties>;

  int32 numMaterial_ = 0;
  ViewType1D<real> rho_;                  // 材料密度数组
  PropertiesArrayType supProperties_; 		// 材料属性对称数组
  real globalScaleFactor_ = 1.0;          // 全局缩放因子

  /**
   * @brief 读取SUP&JKR-F字典参数
   */
  bool readsupDictionary(const dictionary &dict) {
    // 读取全局缩放因子
    if (dict.containsDataEntry("scaleFactor")) {
      globalScaleFactor_ = dict.getVal<real>("scaleFactor");
    } else {
      globalScaleFactor_ = 1.0;
    }

    // 读取材料参数向量
    auto Yeff = dict.getVal<realVector>("Yeff");   	// 有效杨氏模量
    auto Geff = dict.getVal<realVector>("Geff");   	// 有效剪切模量
    auto nu = dict.getVal<realVector>("nu");       	// 泊松比
    auto gamma = dict.getVal<realVector>("gamma"); 	// 表面能密度
    auto en = dict.getVal<realVector>("en");       	// 法向恢复系数
    auto mu = dict.getVal<realVector>("mu");       	// 摩擦系数

    // 可选参数
    realVector et;
    if (dict.containsDataEntry("et")) {
      et = dict.getVal<realVector>("et");
    } else {
      et = en; // 默认切向恢复系数等于法向恢复系数
    }

    auto nElem = Yeff.size();

    // 验证尺寸一致性
    if (nElem != nu.size() || nElem != gamma.size() || nElem != en.size() ||
        nElem != mu.size()) {
      fatalErrorInFunction << "Size mismatch in material properties vectors.\n";
      return false;
    }

    // 检查对称数组尺寸
    uint32 nMat;
    if (!PropertiesArrayType::getN(nElem, nMat)) {
      fatalErrorInFunction
          << "Properties sizes do not match a symmetric array.\n";
      return false;
    }

    if (numMaterial_ != nMat) {
      fatalErrorInFunction << "Expected " << numMaterial_
                           << " materials but got " << nMat << "\n";
      return false;
    }

		// TODO: 以上改动后续需要验证
    // 创建属性数组
    Vector<supProperties> prop("prop", nElem);
    ForAll(i, Yeff) {
      prop[i] = {
          globalScaleFactor_, Yeff[i], Geff[i], gamma[i], en[i], et[i], mu[i]};
    }

    supProperties_.assign(prop);
    return true;
  }

  static const char *modelName() {
    if constexpr (limited)
      return "supLimited";
    else
      return "supNonLimited";
  }

public:
  TypeInfoNV(modelName());

  INLINE_FUNCTION_HD
  sup() {}

  /**
   * @brief 构造函数
   */
  sup(int32 nMaterial, const ViewType1D<real> &rho, const dictionary &dict)
      : numMaterial_(nMaterial), rho_("rho", nMaterial),
        supProperties_("supProperties", nMaterial) {
    Kokkos::deep_copy(rho_, rho);
    if (!readsupDictionary(dict)) {
      fatalExit;
    }
  }

  INLINE_FUNCTION_HD
  sup(const sup &) = default;

  INLINE_FUNCTION_HD
  sup(sup &&) = default;

  INLINE_FUNCTION_HD
  sup &operator=(const sup &) = default;

  INLINE_FUNCTION_HD
  sup &operator=(sup &&) = default;

  INLINE_FUNCTION_HD
  ~sup() = default;

  INLINE_FUNCTION_HD
  int32 numMaterial() const { return numMaterial_; }

  /**
   * @brief SUP&JKR-F接触力计算
   *
   * 算法步骤：
   * 1. 将缩放粒子变量转换到原始粒子尺度
   * 2. 在原始粒子尺度计算JKR-F接触力
   * 3. 将力缩放回缩放粒子系统（l²缩放）
   *
   * JKR-F模型特点：
   * - 使用表面能密度γ而非体积能密度
   * - 拉脱力F_po = 3πγR_eff（与完整JKR模型相同）
   * - 在δ=0时接触激活/失活（无需处理负重叠）
   */
  INLINE_FUNCTION_HD
  void contactForce(const real dt, const uint32 i, const uint32 j,
                    const uint32 propId_i, const uint32 propId_j,
                    const real Ri,      // 缩放粒子半径
                    const real Rj,      // 缩放粒子半径
                    const real ovrlp_n, // 缩放粒子重叠量
                    const realx3 &Vr,   // 相对速度
                    const realx3 &Nij,  // 法向单位向量（从i指向j）
                    contactForceStorage &history,
                    realx3 &FCn, // 输出：法向力（缩放尺度）
                    realx3 &FCt  // 输出：切向力（缩放尺度）
  ) const {
    // 获取材料属性
    auto prop = supProperties_(propId_i, propId_j);
    const real l = prop.l_;         // 缩放因子
    const real gamma = prop.gamma_; // 表面能密度

    // ========== 步骤1：变量转换到原始粒子尺度 ==========
    const real Ri_o = Ri / l;                            // 原始粒子半径
    const real Rj_o = Rj / l;                            // 原始粒子半径
    const real Reff_o = 1.0 / (1.0 / Ri_o + 1.0 / Rj_o); // 有效半径
    const real delta_o = ovrlp_n / l; // 原始粒子重叠量

    // JKR-F模型在δ=0时激活/失活
    if (delta_o <= 0.0) {
      FCn = 0.0;
      FCt = 0.0;
      history.overlap_t_ = 0.0;
      return;
    }

    // 速度分解
    const real vrn = dot(Vr, Nij);    // 法向相对速度
    const realx3 Vt = Vr - vrn * Nij; // 切向相对速度

    // ========== 步骤2：计算原始粒子的JKR-F接触力 ==========

    // 计算质量（原始粒子）
    const real mi_o = 3 * Pi / 4 * pow(Ri_o, 3) * rho_[propId_i]; // 使用近似公式
    const real mj_o = 3 * Pi / 4 * pow(Rj_o, 3) * rho_[propId_j];
    const real meff_o = (mi_o * mj_o) / (mi_o + mj_o);

    // JKR-F法向力（原始粒子尺度）
    // F_n = (4/3)√(R_eff*E*)δ^(3/2) - √(16πγE*)(R_eff*δ)^(3/4)
    const real hertz_term =
        (4.0 / 3.0) * sqrt(Reff_o * prop.Yeff_) * pow(delta_o, 1.5);
		
		// 弹性力初始化为Hertz力
		real Fn_elastic_o = hertz_term;

		// 只在有黏附时计算黏附项
		if (gamma > zero) {
			// 弹性力（压缩为正）- 黏附力（拉伸）	
			const real adhesion_term = 
					sqrt(16.0 * Pi * gamma * prop.Yeff_) * pow(Reff_o * delta_o, 0.75);
				
			Fn_elastic_o -= adhesion_term;  // 减去黏附力
		}

    // 计算法向阻尼
    real Fn_damping_o = 0.0;
    if (prop.en_ < 1.0 && prop.en_ > 0.0 && vrn != 0.0) {
      // Hertz接触刚度用于阻尼计算
      const real K_hertz = (4.0 / 3.0) * prop.Yeff_ * sqrt(Reff_o);
      const real sqrt_meff_K = sqrt(meff_o * K_hertz);

      // Tsuji阻尼系数
      const real ethan =
          -2.2664 * log(prop.en_) / sqrt(pow(log(prop.en_), 2) + pow(Pi, 2));

      Fn_damping_o = sqrt_meff_K * ethan * pow(delta_o, 0.25) * vrn;
    }

    // 总法向力（原始粒子尺度）
    realx3 FCn_o = (Fn_elastic_o - Fn_damping_o) * Nij;

    // ========== 切向力计算（原始粒子尺度）==========

    // 更新切向重叠历史（从缩放尺度转换）
    realx3 overlap_t_o = history.overlap_t_ / l;
    overlap_t_o += Vt * dt;

    // 切向刚度（Mindlin理论）
    const real kt_hertz = 8.0 * prop.Geff_ * sqrt(Reff_o * delta_o);

    // 切向弹性力
    realx3 FCt_o = -kt_hertz * overlap_t_o;

    // 切向阻尼（可选）
    if (prop.et_ < 1.0 && prop.et_ > 0.0 && length(Vt) > zero) {
      real ethat = -2.0 * log(prop.et_) * sqrt(kt_hertz * meff_o) /
                   sqrt(pow(log(prop.et_), 2) + pow(Pi, 2));
      FCt_o -= ethat * Vt;
    }

    // ========== 摩擦力限制（考虑黏附的Coulomb定律）==========

		// 默认摩擦极限（无黏附）
		const real fn_mag = length(FCn_o);
		real ft_limit = prop.mu_ * fn_mag;

		// 只在有黏附时修正摩擦极限
		if (gamma > zero) {
			// JKR-F拉脱力（原始粒子尺度）
			const real Fpo_o = 3.0 * Pi * gamma * Reff_o;
			// 使用Thornton方法：摩擦极限基于(|Fn| + 2*Fpo)
			ft_limit = prop.mu_ * (fn_mag + 2.0 * Fpo_o);
		}

		const real ft_mag = length(FCt_o);

    if (ft_mag > ft_limit) {
      if (length(overlap_t_o) > zero) {
        if constexpr (limited) {
          // 限制模式：调整力并更新切向重叠
          FCt_o *= (ft_limit / ft_mag);
          overlap_t_o = -(FCt_o / kt_hertz);
        } else {
          // 非限制模式：仅调整力
          FCt_o = (FCt_o / ft_mag) * ft_limit;
        }
      } else {
        FCt_o = 0.0;
        overlap_t_o = 0.0;
      }
    }

    // 更新历史（转换回缩放粒子坐标）
    history.overlap_t_ = overlap_t_o * l;

    // ========== 步骤3：力缩放到缩放粒子系统 ==========
    // SUP模型：F_S = l² * F_O
    const real l_squared = l * l;
    FCn = FCn_o * l_squared;
    FCt = FCt_o * l_squared;
  }

  /**
   * @brief 获取拉脱力（缩放粒子尺度）
   *
   * JKR-F模型的拉脱力：F_po = 3πγR_eff
   * 缩放后：F_po_S = l² * F_po_O
   */
  INLINE_FUNCTION_HD
  real getPullOffForce(const uint32 propId_i, const uint32 propId_j,
                       const real Ri, // 缩放粒子半径
                       const real Rj  // 缩放粒子半径
  ) const {
    auto prop = supProperties_(propId_i, propId_j);
    const real l = prop.l_;
    const real gamma = prop.gamma_;

		if (gamma <= zero) {
			return 0.0;
		}

    // 原始粒子有效半径
    const real Ri_o = Ri / l;
    const real Rj_o = Rj / l;
    const real Reff_o = 1.0 / (1.0 / Ri_o + 1.0 / Rj_o);

    // 原始粒子拉脱力
    const real Fpo_o = 3.0 * Pi * gamma * Reff_o;

    // 缩放到缩放粒子系统
    return Fpo_o * l * l;
  }
};

/**
 * @brief 带CDT滚动阻力的SUP&JKR-F模型
 */
template <bool limited = true> class sup_CDT : public sup<limited> {
public:
  using contactForceStorage = typename sup<limited>::contactForceStorage;

protected:
  realSymArray_D mur_; // 滚动摩擦系数矩阵

  bool readRollingDict(const dictionary &dict) {
    auto mur = dict.getVal<realVector>("mur");

    uint32 nMat;
    if (!realSymArray_D::getN(mur.size(), nMat) ||
        nMat != this->numMaterial()) {
      fatalErrorInFunction << "Wrong number of values in mur.\n";
      return false;
    }

    mur_.assign(mur);
    return true;
  }

public:
  TypeInfoNV(word("sup_CDT<" + sup<limited>::TYPENAME() + ">"));

  sup_CDT(int32 nMaterial, const ViewType1D<real> &rho,
              const dictionary &dict)
      : sup<limited>(nMaterial, rho, dict), mur_("mur", nMaterial) {
    if (!readRollingDict(dict)) {
      fatalExit;
    }
  }

  /**
   * @brief 计算滚动阻力矩（考虑黏附）
   */
  INLINE_FUNCTION_HD
  void rollingFriction(const real dt, const uint32 i, const uint32 j,
                       const uint32 propId_i, const uint32 propId_j,
                       const real Ri, const real Rj, const realx3 &wi,
                       const realx3 &wj, const realx3 &Nij,
                       const realx3 &FCn, // 法向力（缩放尺度）
                       realx3 &Mri, realx3 &Mrj) const {
    // 获取缩放参数
    auto prop = this->supProperties_(propId_i, propId_j);
    const real l = prop.l_;
    const real gamma = prop.gamma_;

    // 原始粒子半径
    const real Ri_o = Ri / l;
    const real Rj_o = Rj / l;
    const real Reff_o = (Ri_o * Rj_o) / (Ri_o + Rj_o);

    // 角速度转换：ω_O = l × ω_S
    const realx3 wi_o = wi * l;
    const realx3 wj_o = wj * l;
    const realx3 w_rel_o = wi_o - wj_o;
    const real w_mag = length(w_rel_o);

    if (w_mag > 1e-10) {
			// 计算法向力大小（原始尺度）
			const real Fn_s_mag = length(FCn);
			const real Fn_o_mag = Fn_s_mag / (l * l);
			
			// 有效法向力（用于滚动阻力）
			real Fn_eff_o = Fn_o_mag;
			
			// 只在有黏附时添加拉脱力贡献
			if (gamma > zero) {
					const real Fpo_o = 3.0 * Pi * gamma * Reff_o;
					Fn_eff_o = Fn_o_mag + 2.0 * Fpo_o;
			}
			
			// CDT滚动阻力矩（原始尺度）
			realx3 Mr_o = -mur_(propId_i, propId_j) * Reff_o * Fn_eff_o * (w_rel_o / w_mag);
			
			// 缩放到缩放粒子系统
			Mri = Mr_o * (l * l);
			Mrj = -Mri;
		} else {
      Mri = 0.0;
      Mrj = 0.0;
    }
  }
};

} // namespace pFlow::cfModels

#endif // __supCF_hpp__