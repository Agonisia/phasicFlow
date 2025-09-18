#ifndef __supRollingCDT_hpp__
#define __supRollingCDT_hpp__

namespace pFlow::cfModels
{

/**
 * @brief SUP模型的CDT滚动阻力实现
 * 
 * CDT (Constant Directional Torque) 模型：滚动阻力矩的方向始终与相对角速度相反
 * 
 * SUP缩放关系：
 * - 角速度: ω_O = l × ω_S (原始粒子角速度 = 缩放因子 × 缩放粒子角速度)
 * - 力矩: M_S = l² × M_O (缩放粒子力矩 = 缩放因子² × 原始粒子力矩)
 * 
 * @tparam contactForceModel 基础接触力模型（应该是sup<true>或sup<false>）
 */
template<typename contactForceModel>
class supRollingCDT : public contactForceModel
{
public:
	using contactForceStorage = typename contactForceModel::contactForceStorage;
    
// protected:
	realSymArray_D mur_;  // 滚动摩擦系数矩阵

	/**
	 * @brief 读取滚动阻力相关参数
	 */
	bool readSupRollingDict(const dictionary& dict)
	{
		auto mur = dict.getVal<realVector>("mur");
		
		uint32 nMat;
		if(!realSymArray_D::getN(mur.size(), nMat) || nMat != this->numMaterial())
		{
			fatalErrorInFunction<<
			"wrong number of values supplied in mur.\n";
			return false;
		}
		
		mur_.assign(mur);
		return true;
	}

public:
	// 类型信息，用于运行时识别
	TypeInfoNV(word("supRollingCDT<"+contactForceModel::TYPENAME()+">"));
	

	/**
	 * @brief 构造函数
	 */
	supRollingCDT(int32 nMaterial, const ViewType1D<real>& rho, const dictionary& dict)
	:
		contactForceModel(nMaterial, rho, dict),
		mur_("mur", nMaterial)
	{
		if(!readSupRollingDict(dict))
		{
			fatalExit;
		}
	}
	
	/**
	 * @brief 计算滚动阻力矩
	 * 
	 * 算法步骤：
	 * 1. 将缩放粒子的角速度转换为原始粒子角速度
	 * 2. 在原始粒子尺度下计算滚动阻力矩
	 * 3. 将力矩缩放回缩放粒子系统（l²缩放）
	 * 
	 * @param dt 时间步长
	 * @param i, j 粒子索引
	 * @param propId_i, propId_j 材料属性索引
	 * @param Ri, Rj 缩放粒子半径
	 * @param wi, wj 缩放粒子角速度
	 * @param Nij 接触法向单位向量
	 * @param FCn 法向力（缩放尺度，包含弹性和粘附分量）
	 * @param Mri, Mrj 输出的滚动阻力矩（缩放尺度）
	 */
	INLINE_FUNCTION_HD
	void rollingFriction(
		const real dt,
		const uint32 i,
		const uint32 j,
		const uint32 propId_i,
		const uint32 propId_j,
		const real Ri,          // 缩放粒子半径
		const real Rj,          // 缩放粒子半径
		const realx3& wi,       // 缩放粒子角速度
		const realx3& wj,       // 缩放粒子角速度
		const realx3& Nij,      // 接触法向单位向量
		const realx3& FCn,      // 法向力（缩放尺度）
		realx3& Mri,            // 输出：粒子i的滚动阻力矩
		realx3& Mrj             // 输出：粒子j的滚动阻力矩
	)const
	{
		// ========== 步骤1：获取SUP缩放参数 ==========
		auto prop = this->supProperties_(propId_i, propId_j);
		const real l = prop.l_;  // 缩放因子 l = d_S/d_O
		
		// ========== 步骤2：转换到原始粒子尺度 ==========
		// 原始粒子半径
		const real Ri_o = Ri / l;
		const real Rj_o = Rj / l;
		
		// 原始粒子有效半径
		const real Reff_o = (Ri_o * Rj_o) / (Ri_o + Rj_o);
		
		// 角速度转换：ω_O = l × ω_S
		// 注意：这是SUP模型的关键转换关系
		const realx3 wi_o = wi * l;  // 原始粒子i的角速度
		const realx3 wj_o = wj * l;  // 原始粒子j的角速度
		
		// 相对角速度（原始尺度）
		const realx3 w_rel_o = wi_o - wj_o;
		const real w_mag = length(w_rel_o);
		
		// ========== 步骤3：计算原始粒子的滚动阻力矩 ==========
		if(w_mag > 1e-10)  // 避免除零
		{
			// 法向力大小（原始尺度）
			// 注意：FCn是缩放尺度的力，需要转换
			// F_S = l² × F_O，所以 F_O = F_S / l²
			const real Fn_o_mag = length(FCn) / (l * l);
			
			// CDT模型：滚动阻力矩与相对角速度方向相反
			// M_R = -μ_r × R_eff × F_N × (ω_rel/|ω_rel|)
			realx3 Mr_o = -mur_(propId_i, propId_j) * Reff_o * Fn_o_mag * 
										(w_rel_o / w_mag);
			
			// ========== 步骤4：力矩缩放回缩放粒子系统 ==========
			// SUP缩放：M_S = l² × M_O
			Mri = Mr_o * (l * l);
			
			// 作用-反作用原理
			Mrj = -Mri;
		}
		else
		{
			// 无相对角速度，无滚动阻力
			Mri = 0.0;
			Mrj = 0.0;
		}
		
		// 注意：这里不需要移除法向分量，因为CDT模型中
		// 滚动阻力矩的方向是相对角速度的方向，已经自然地在接触平面内
	}
};

} // namespace pFlow::cfModels

#endif // __supRollingCDT_hpp__