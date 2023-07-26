for i in 1:5
    t = @elapsed begin
    sleep(0.3)
    end 
    print(t)
end