function [Q,R] = qr_own(A)
    % Vlastní QR rozklad pomocí Gramova--Schmidtova procesu.
    % Předpokládá lineárně nezávislé sloupce matice A.

    [~,n] = size(A);       % Počet sloupců matice A
    R = eye(n);            % Koeficienty rozkladu budou v horní
                           % trojúhelníkové matici R
    Q = A;                 % Sloupce Q postupně ortogonalizujeme

    % Ortogonalizace sloupců matice A
    for j = 2:n
        for i = 1:j-1
            % Koeficient projekce j-tého sloupce na i-tý sloupec Q.
            % Sloupce Q jsou zatím ortogonální, ale nenormalizované.
            R(i,j) = dot(Q(:,i),A(:,j)) / dot(Q(:,i),Q(:,i));

            % Odečtením projekce odstraníme složku ve směru Q(:,i).
            Q(:,j) = Q(:,j) - R(i,j)*Q(:,i);
        end
    end

    % Normy dosud ortogonálních, ale nenormalizovaných sloupců Q
    d = sqrt(diag(Q'*Q));

    % Normalizujeme sloupce Q, aby platilo Q'*Q = I
    Q = Q*diag(1./d);

    % Upravíme R tak, aby i po normalizaci stále platilo A = Q*R
    R = diag(d)*R;
end