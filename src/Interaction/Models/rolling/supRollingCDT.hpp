#ifndef __supRollingCDT_hpp__
#define __supRollingCDT_hpp__

namespace pFlow::cfModels
{

template<typename contactForceModel>
class supRollingCDT : public contactForceModel
{
public:
    using contactForceStorage = typename contactForceModel::contactForceStorage;
    
    realSymArray_D mur_;
    real scaleFactor_;

public:
    TypeInfoNV(word("supRollingCDT<"+contactForceModel::TYPENAME()+">"));
    
    supRollingCDT(int32 nMaterial, const ViewType1D<real>& rho, const dictionary& dict)
    :
        contactForceModel(nMaterial, rho, dict),
        mur_("mur", nMaterial)
    {
        auto mur = dict.getVal<realVector>("mur");
        mur_.assign(mur);
        scaleFactor_ = dict.getVal<real>("scaleFactor");
    }
    
    INLINE_FUNCTION_HD
    void rollingFriction(
        const real dt,
        const uint32 i,
        const uint32 j,
        const uint32 propId_i,
        const uint32 propId_j,
        const real Ri,
        const real Rj,
        const realx3& wi,
        const realx3& wj,
        const realx3& Nij,
        const realx3& FCn,
        realx3& Mri,
        realx3& Mrj
    )const
    {
        // CDT模型: M_R = -μ_r * r* * F_N * (ω_rel / |ω_rel|)
        
        // 原始粒子参数
        real l = scaleFactor_;
        real Ri_o = Ri / l;
        real Rj_o = Rj / l;
        real Reff_o = (Ri_o * Rj_o) / (Ri_o + Rj_o);
        
        // 原始粒子相对角速度 (ω_O = l×ω_S)
        realx3 w_rel_o = (wi - wj) * l;
        real w_mag = length(w_rel_o);
        
        if(w_mag > 1e-10)
        {
            // 原始粒子法向力
            real Fn_o = length(FCn) / (l * l);
            
            // CDT滚动力矩（原始尺度）
            realx3 Mr_o = -mur_(propId_i, propId_j) * Reff_o * Fn_o * 
                          (w_rel_o / w_mag);
            
            // 缩放回缩放粒子系统 (M_S = l²×M_O)
            Mri = Mr_o * (l * l);
            Mrj = -Mri;
        }
        else
        {
            Mri = Mrj = realx3(0.0);
        }
    }
};

} // namespace pFlow::cfModels

#endif // __supRollingCDT_hpp__