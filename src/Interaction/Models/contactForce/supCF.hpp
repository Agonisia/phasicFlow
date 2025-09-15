// src/Interaction/Models/contactForce/supCF.hpp
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
        realx3 overlap_t_ = 0.0;
    };
    
    struct supProperties
    {
        real scaleFactor_ = 1.0;
        real DE_ = 87e-6;       // 特征长度尺度
        real muEffInf_ = 0.447; // 有效摩擦系数
        real beta_ = 0.014;     // 经验参数
        
        // 材料参数（原始粒子）
        real Yeff_ = 6.3e10;
        real Geff_;
        real surfaceEnergy_ = 0.05;
        real en_ = 0.8;
        real mu_ = 0.5;
        
        INLINE_FUNCTION_HD
        supProperties() {
            // 计算Geff从Yeff和泊松比
            real nu = 0.24;
            Geff_ = Yeff_ / (2.0 * (1.0 + nu));
        }
        
        INLINE_FUNCTION_HD
        supProperties(real scaleFactor, real DE, real Yeff, real en, real mu):
            scaleFactor_(scaleFactor), DE_(DE), Yeff_(Yeff), en_(en), mu_(mu) 
        {
            real nu = 0.24;
            Geff_ = Yeff_ / (2.0 * (1.0 + nu));
        }
    };

protected:
    using SUPArrayType = symArray<supProperties>;
    
    int32 numMaterial_ = 0;
    ViewType1D<real> rho_;
    SUPArrayType supProperties_;
    
    bool readSUPDictionary(const dictionary& dict)
    {
        // 读取SUP模型参数
        real scaleFactor = dict.getVal<real>("scaleFactor");
        real DE = dict.getValOrDefault<real>("DE", 87e-6);
        
        auto Yeff = dict.getVal<realVector>("Yeff");
        auto en = dict.getVal<realVector>("en");
        auto mu = dict.getVal<realVector>("mu");
        
        auto nElem = Yeff.size();
        
        // 验证尺寸
        uint32 nMat;
        if(!SUPArrayType::getN(nElem, nMat) || nMat != numMaterial_)
        {
            fatalErrorInFunction
            "Size mismatch for SUP properties.\n";
            return false;
        }
        
        Vector<supProperties> prop("prop", nElem);
        ForAll(i, Yeff)
        {
            prop[i] = {scaleFactor, DE, Yeff[i], en[i], mu[i]};
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
    int32 numMaterial()const { return numMaterial_; }
    
    // 核心SUP模型计算
    INLINE_FUNCTION_HD
    void contactForce(
        const real dt,
        const uint32 i,
        const uint32 j,
        const uint32 propId_i,
        const uint32 propId_j,
        const real Ri,
        const real Rj,
        const real ovrlp_n,
        const realx3& Vr,
        const realx3& Nij,
        contactForceStorage& history,
        realx3& FCn,
        realx3& FCt
    )const
    {
        auto prop = supProperties_(propId_i, propId_j);
        real l = prop.scaleFactor_;
        
        // 1. 转换到原始粒子变量
        real R_o = (Ri + Rj) / (2.0 * l);  // 原始半径
        real delta_o = ovrlp_n / l;         // 原始重叠量
        realx3 V_o = Vr;                    // v_O = v_S (速度不变)
        
        // 2. 计算原始粒子的有效参数
        real Reff_o = R_o / 2.0;  // 对于相同大小粒子
        
        // 3. 计算原始粒子的力（使用Hertz-JKR模型）
        real K_hertz = 4.0/3.0 * prop.Yeff_ * sqrt(Reff_o);
        
        // 法向力（原始尺度）
        real vrn = dot(V_o, Nij);
        
        // 质量（原始粒子）
        real m_o = 4.0/3.0 * Pi * pow(R_o, 3) * rho_[propId_i];
        real meff = m_o / 2.0;  // 对于相同质量粒子
        
        // 阻尼系数
        real ethan = -2.0 * log(prop.en_) * sqrt(K_hertz * meff) /
                     sqrt(pow(log(prop.en_), 2) + pow(Pi, 2));
        
        // 原始粒子的法向力
        realx3 FCn_o = (-K_hertz * pow(delta_o, 1.5) - 
                        ethan * pow(delta_o, 0.25) * vrn) * Nij;
        
        // Van der Waals力（原始尺度）
        real A_H = 24 * Pi * 1.65e-20 * prop.surfaceEnergy_;  // Hamaker常数
        realx3 F_vdw_o = (A_H * Reff_o / (6.0 * pow(1.65e-10, 2))) * Nij;
        
        FCn_o = FCn_o + F_vdw_o;
        
        // 切向力（原始尺度）
        realx3 Vt = Vr - vrn * Nij;
        history.overlap_t_ += Vt * dt / l;  // 调整切向重叠量
        
        real kt_o = 8.0 * prop.Geff_ * sqrt(Reff_o * delta_o);
        realx3 FCt_o = -kt_o * history.overlap_t_;
        
        // 摩擦力限制
        real ft = length(FCt_o);
        real ft_fric = prop.mu_ * length(FCn_o);
        
        if(ft > ft_fric)
        {
            if(length(history.overlap_t_) > 0.0)
            {
                if constexpr (limited)
                {
                    FCt_o *= (ft_fric/ft);
                    history.overlap_t_ = -FCt_o/kt_o;
                }
                else
                {
                    FCt_o = (FCt_o/ft) * ft_fric;
                }
            }
            else
            {
                FCt_o = 0.0;
            }
        }
        
        // 4. 缩放力回到缩放粒子系统（l²缩放）
        FCn = FCn_o * (l * l);
        FCt = FCt_o * (l * l);
    }
};

} // namespace pFlow::cfModels

#endif // __supCF_hpp__