 function v=shuffle(v)
    if isrow(v)
     v=v(randperm(length(v)));
    else
      v=v(randperm(size(v,1)),:);
    end
 end