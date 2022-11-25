macro m(ex)
    protected_ex = QuoteNode(ex)
    :(println($protected_ex, " = ", $(esc(ex))))
end