#ifndef __supCF_hpp__
#define __supCF_hpp__

#include "types.hpp"
#include "symArrays.hpp"

namespace pFlow::cfModels
{

template<bool limited=true>
class sup
{
public:
    struct contactForceStorage
    {
        realx3 overlap_t_ = 0.0;  // 切向重叠历史（缩放粒子坐标系）
    };
    
    struct supProperties
    {
        // 缩放因子
        real l_ = 1.0;           // l = d_S/d_O
        
        // 材料参数（每种材料对）
        real Yeff_ = 1e7;        // 有效杨氏模量 [Pa]
        real Geff_ = 4e6;        // 有效剪切模量 [Pa]
        real en_ = 0.8;          // 法向恢复系数
        real et_ = 0.8;          // 切向恢复系数（可选）
        real mu_ = 0.5;          // 滑动摩擦系数
        
        INLINE_FUNCTION_HD
        supProperties(){}
        
        INLINE_FUNCTION_HD
        supProperties(real l, real Yeff, real Geff, real en, real mu):
            l_(l), Yeff_(Yeff), Geff_(Geff), en_(en), et_(en), mu_(mu) {}
        
        INLINE_FUNCTION_HD
        supProperties(const supProperties&) = default;
        
        INLINE_FUNCTION_HD
        supProperties& operator=(const supProperties&) = default;
        
        INLINE_FUNCTION_HD
        ~supProperties() = default;
    };

protected:
    using SUPArrayType = symArray<supProperties>;
    
    int32 numMaterial_ = 0;
    ViewType1D<real> rho_;
    SUPArrayType supProperties_;
    real globalScaleFactor_ = 1.0;
    
    bool readSUPDictionary(const dictionary& dict)
    {
        // 读取全局缩放因子
        globalScaleFactor_ = dict.getValOrDefault<real>("scaleFactor", 1.0);
        
        // 读取材料参数向量
        auto Yeff = dict.getVal<realVector>("Yeff");
        auto Geff = dict.getVal<realVector>("Geff");  
        auto nu = dict.getVal<realVector>("nu");       // 泊松比
        auto en = dict.getVal<realVector>("en");
        auto mu = dict.getVal<realVector>("mu");
        
        auto nElem = Yeff.size();
        
        // 验证尺寸一致性
        if(nElem != nu.size())
        {
            fatalErrorInFunction
            "sizes of Yeff("<<nElem<<") and nu("<<nu.size()<<") do not match.\n";
            return false;
        }
        
        if(nElem != en.size())
        {
            fatalErrorInFunction
            "sizes of Yeff("<<nElem<<") and en("<<en.size()<<") do not match.\n";
            return false;
        }
        
        if(nElem != mu.size())
        {
            fatalErrorInFunction
            "sizes of Yeff("<<nElem<<") and mu("<<mu.size()<<") do not match.\n";
            return false;
        }
        
        // 检查对称数组尺寸
        uint32 nMat;
        if(!SUPArrayType::getN(nElem, nMat))
        {
            fatalErrorInFunction
            "sizes of properties do not match a symmetric array.\n";
            return false;
        }
        
        if(numMaterial_ != nMat)
        {
            fatalErrorInFunction
            "size mismatch for properties. Expected "<<numMaterial_
            " materials but got "<<nMat<<"\n";
            return false;
        }
        
        // 如果Geff未提供，从Yeff和nu计算
        if(Geff.size() != nElem)
        {
            Geff.resize(nElem);
            ForAll(i, Yeff)
            {
                // G = E / (2(1+ν))
                Geff[i] = Yeff[i] / (2.0 * (1.0 + nu[i]));
            }
        }
        
        // 创建属性数组
        Vector<supProperties> prop("prop", nElem);
        ForAll(i, Yeff)
        {
            prop[i] = {globalScaleFactor_, Yeff[i], Geff[i], en[i], mu[i]};
        }
        
        supProperties_.assign(prop);
        return true;
    }
    
    static const char* modelName()
    {
        if constexpr (limited)
            return "supLimited";
        else
            return "supNonLimited";
    }

public:
    TypeInfoNV(modelName());
    
    INLINE_FUNCTION_HD
    sup(){}
    
    sup(int32 nMaterial, const ViewType1D<real>& rho, const dictionary& dict)
    :
        numMaterial_(nMaterial),
        rho_("rho", nMaterial),
        supProperties_("supProperties", nMaterial)
    {
        Kokkos::deep_copy(rho_, rho);
        if(!readSUPDictionary(dict))
        {
            fatalExit;
        }
    }
    
    INLINE_FUNCTION_HD
    sup(const sup&) = default;
    
    INLINE_FUNCTION_HD
    sup(sup&&) = default;
    
    INLINE_FUNCTION_HD
    sup& operator=(const sup&) = default;
    
    INLINE_FUNCTION_HD
    sup& operator=(sup&&) = default;
    
    INLINE_FUNCTION_HD
    ~sup() = default;
    
    INLINE_FUNCTION_HD
    int32 numMaterial()const { return numMaterial_; }
    
    /**
     * SUP-Hertz接触力计算
     * 遵循SUP模型的三步流程：
     * 1. 缩放粒子变量 -> 原始粒子变量
     * 2. 在原始粒子尺度计算Hertz接触力
     * 3. 力缩放回缩放粒子系统（l²缩放）
     */
    INLINE_FUNCTION_HD
    void contactForce(
        const real dt,
        const uint32 i,
        const uint32 j,
        const uint32 propId_i,
        const uint32 propId_j,
        const real Ri,          // 缩放粒子半径
        const real Rj,          // 缩放粒子半径
        const real ovrlp_n,     // 缩放粒子重叠量
        const realx3& Vr,       // 相对速度（v_S = v_O）
        const realx3& Nij,      // 法向单位向量
        contactForceStorage& history,
        realx3& FCn,            // 输出：法向力（缩放尺度）
        realx3& FCt             // 输出：切向力（缩放尺度）
    )const
    {
        // 获取材料属性
        auto prop = supProperties_(propId_i, propId_j);
        const real l = prop.l_;  // 缩放因子
        
        // ========== 步骤1：变量转换到原始粒子尺度 ==========
        const real Ri_o = Ri / l;                        // 原始粒子i半径
        const real Rj_o = Rj / l;                        // 原始粒子j半径
        const real Reff_o = 1.0 / (1.0/Ri_o + 1.0/Rj_o); // 有效半径
        const real delta_o = ovrlp_n / l;                // 原始粒子重叠量
        
        // 检查是否有接触
        if(delta_o <= 0.0)
        {
            FCn = 0.0;
            FCt = 0.0;
            history.overlap_t_ = 0.0;
            return;
        }
        
        // 速度分解
        const real vrn = dot(Vr, Nij);      // 法向相对速度
        const realx3 Vt = Vr - vrn * Nij;   // 切向相对速度
        
        // ========== 步骤2：计算原始粒子的Hertz接触力 ==========
        
        // 计算质量和有效质量（原始粒子）
        const real mi_o = (4.0/3.0) * Pi * pow(Ri_o, 3) * rho_[propId_i]; // 注意这里使用的更精确的体积公式
        const real mj_o = (4.0/3.0) * Pi * pow(Rj_o, 3) * rho_[propId_j];
        const real meff_o = (mi_o * mj_o) / (mi_o + mj_o);
        
        // Hertz接触刚度
        const real K_hertz = (4.0/3.0) * prop.Yeff_ * sqrt(Reff_o);
        const real sqrt_meff_K = sqrt(meff_o * K_hertz);
        
        // 法向阻尼系数（基于恢复系数）
        real ethan = 0.0;
        if(prop.en_ < 1.0 && prop.en_ > 0.0)
        {
            // Tsuji等人的阻尼模型
            ethan = -2.2664 * log(prop.en_) / 
                    sqrt(pow(log(prop.en_), 2) + pow(Pi, 2));
        }
        
        // Hertz法向力（原始粒子尺度）
        realx3 FCn_o = (
            -K_hertz * pow(delta_o, 1.5) -                    // 弹性力
            sqrt_meff_K * ethan * pow(delta_o, 0.25) * vrn    // 阻尼力
        ) * Nij;
        
        // ========== 切向力计算（原始粒子尺度）==========
        
        // 更新切向重叠历史
        // 注意：历史以缩放粒子坐标存储，需要转换
        realx3 overlap_t_o = history.overlap_t_ / l;  // 转换到原始尺度
        overlap_t_o += Vt * dt;                        // 累积切向位移
        
        // 切向刚度（Mindlin理论）
        const real kt_hertz = 8.0 * prop.Geff_ * sqrt(Reff_o * delta_o);
        
        // 切向弹性力
        realx3 FCt_o = -kt_hertz * overlap_t_o;
        
        // 切向阻尼（可选）
        if(prop.et_ < 1.0 && prop.et_ > 0.0)
        {
            real ethat = -2.0 * log(prop.et_) * sqrt(kt_hertz * meff_o) /
                        sqrt(pow(log(prop.et_), 2) + pow(Pi, 2));
            FCt_o -= ethat * Vt;
        }
        
        // ========== 摩擦力限制（Coulomb定律）==========
        const real ft = length(FCt_o);
        const real ft_fric = prop.mu_ * length(FCn_o);
        
        if(ft > ft_fric)
        {
            if(length(overlap_t_o) > zero)
            {
                if constexpr (limited)
                {
                    // 限制模式：调整力并更新切向重叠
                    FCt_o *= (ft_fric / ft);
                    overlap_t_o = -(FCt_o / kt_hertz);
                }
                else
                {
                    // 非限制模式：仅调整力
                    FCt_o = (FCt_o / ft) * ft_fric;
                }
            }
            else
            {
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
};

} // namespace pFlow::cfModels

#endif // __supCF_hpp__