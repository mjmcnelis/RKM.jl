
function make_code_name(name)
    name_split = split(string(name), "_")
    code_name = "" 
    for person in filter(x -> x != "", filter.(!isdigit, name_split))
        code_name *= person[1]
    end
    for order in filter(x -> x != "", filter.(isdigit, name_split))
        code_name *= order
    end
    code_name
end