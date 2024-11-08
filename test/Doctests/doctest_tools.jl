
module MDDocTestTools

using Markdown

export extract_section

function rerender_markdown_subsection(md::Markdown.MD, start_idx, end_idx)
    Markdown.plain(md.content[start_idx:end_idx])
    #buffer = IOBuffer()
    #println(md.content[start_idx:end_idx])
    #Markdown.plain(buffer,md.content[start_idx:end_idx])
    #res = String(take!(buffer))
    #println(res)
    #res

end

"""
This function extracts the text in a Markdown file
that is just after a header and before the next header.
"""
function extract_section(markdown::Markdown.MD, start_header_text::AbstractString)
    content = markdown.content

    function find_header_idx()
        function are_strings_equal(str1::AbstractString, str2::AbstractString)
            lowercase(str1) == lowercase(str2)
        end

        for idx in eachindex(content)
            element = content[idx]
            if isa(element, Markdown.Header) &&
               are_strings_equal(first(element.text), start_header_text)
                return idx
            end
        end
    end

    function find_next_header_idx(start_idx)
        for idx = (start_idx+1):length(content)
            element = content[idx]
            if isa(element, Markdown.Header)
                return idx
            end
        end
    end


    header_idx = find_header_idx()
    if isnothing(header_idx)
        return
    end

    next_header_idx = find_next_header_idx(header_idx)
    end_idx = if isnothing(next_header_idx)
        length(content)
    else
        next_header_idx - 1
    end

    function concat_text(start_idx, end_idx)
        from = start_idx + 1
        until = end_idx

        rerender_markdown_subsection(markdown, from, until)
    end
    concat_text(header_idx, end_idx)
end



end
