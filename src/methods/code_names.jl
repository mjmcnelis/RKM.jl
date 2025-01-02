
function make_code_name(name)
    name_split = split(string(name), "_")
    code_name = ""
    for person in filter(x -> x != "", filter.(!isdigit, name_split))
        code_name *= person[1]
    end
    for order in filter(x -> x != "", filter.(isdigit, name_split))
        code_name *= order
    end
    return code_name
end

function adaptive_code_label(::AdaptiveStepSize)
    return ""
end

function adaptive_code_label(::Doubling)
    return "D"
end

function adaptive_code_label(::CentralDiff)
    return "M"
end
