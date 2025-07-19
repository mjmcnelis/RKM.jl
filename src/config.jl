
struct RKMConfig{OWS, UC, JM, RFM, EMM, SM, OWP, I}
    ode_wrap_y!::OWS
    update_cache::UC
    state_jacobian::JM
    root_finder::RFM
    eigenmax::EMM
    sensitivity::SM
    ode_wrap_p!::OWP
    interpolator::I
end