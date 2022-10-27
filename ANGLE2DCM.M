function dcm = angle2dcm( r1, r2, r3 )
angles = [r1(:) r2(:) r3(:)];

dcm = zeros(3,3,size(angles,1));
cang = cos( angles );
sang = sin( angles );


dcm(1,1,:) = cang(:,2).*cang(:,1);
dcm(1,2,:) = cang(:,2).*sang(:,1);
dcm(1,3,:) = -sang(:,2);
dcm(2,1,:) = sang(:,3).*sang(:,2).*cang(:,1) - cang(:,3).*sang(:,1);
dcm(2,2,:) = sang(:,3).*sang(:,2).*sang(:,1) + cang(:,3).*cang(:,1);
dcm(2,3,:) = sang(:,3).*cang(:,2);
dcm(3,1,:) = cang(:,3).*sang(:,2).*cang(:,1) + sang(:,3).*sang(:,1);
dcm(3,2,:) = cang(:,3).*sang(:,2).*sang(:,1) - sang(:,3).*cang(:,1);
dcm(3,3,:) = cang(:,3).*cang(:,2);
end