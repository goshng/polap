local meta_values = {}

function Meta(meta)
	io.stderr:write("Meta keys:\n")
	for k, v in pairs(meta) do
		-- io.stderr:write(k .. " = " .. pandoc.utils.stringify(v) .. "\n")
		if type(v) == "table" and v.t == "MetaString" then
			meta_values[k] = v.text
		elseif type(v) == "boolean" then
			meta_values[k] = tostring(v)
		elseif type(v) == "number" then
			meta_values[k] = tostring(v)
		elseif type(v) == "string" then
			meta_values[k] = v
		else
			-- fallback for other types (MetaInlines, MetaBlocks, etc.)
			meta_values[k] = pandoc.utils.stringify(v)
		end
	end
	return meta
end

function replace_meta_variables(str)
	return str:gsub("{{(.-)}}", function(key)
		return meta_values[key] or "{{" .. key .. "}}"
	end)
end

function Str(el)
	-- io.stderr:write("Processing Str: " .. el.text .. "\n")
	el.text = replace_meta_variables(el.text)
	return el
end

function CodeBlock(el)
	el.text = replace_meta_variables(el.text)
	return el
end

function Link(el)
	el.target = replace_meta_variables(el.target)
	return el
end
