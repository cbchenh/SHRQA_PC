function PlotCell(Ub,Lb)
    if size(Ub) ~= size(Lb)
       error('The dimensions of UB and LB are different.') 
    end
    [r, d] = size(Ub);
    if d == 2
        for i = 1:r
            rect = [Lb(i,:) Ub(i,:)-Lb(i,:)];
            rectangle('Position',rect,'LineWidth',3);
            hold on;
        end
    end
    if d >= 3
        if d > 3
           disp('------------Only show the first three dimensions-----------') 
        end
        for i = 1:r
            p1 = [Lb(i,1), Lb(i,2), Lb(i,3)];
            p2 = [Ub(i,1), Lb(i,2), Lb(i,3)];
            p3 = [Ub(i,1), Lb(i,2), Ub(i,3)];
            p4 = [Lb(i,1), Lb(i,2), Ub(i,3)];
            p5 = [Lb(i,1), Ub(i,2), Ub(i,3)];
            p6 = [Ub(i,1), Ub(i,2), Ub(i,3)];
            p7 = [Ub(i,1), Ub(i,2), Lb(i,3)];
            p8 = [Lb(i,1), Ub(i,2), Lb(i,3)];
            plot3([p1(1) p2(1) p3(1) p4(1) p1(1) p8(1) p5(1) p4(1) p3(1) p6(1) p5(1) p8(1) p7(1) p6(1) p3(1) p2(1) p7(1)],[p1(2) p2(2) p3(2) p4(2) p1(2) p8(2) p5(2) p4(2) p3(2) p6(2) p5(2) p8(2) p7(2) p6(2) p3(2) p2(2) p7(2)],[p1(3) p2(3) p3(3) p4(3) p1(3) p8(3) p5(3) p4(3) p3(3) p6(3) p5(3) p8(3) p7(3) p6(3) p3(3) p2(3) p7(3)],'LineWidth',3);
            hold on;
        end
    end
    hold off;
end