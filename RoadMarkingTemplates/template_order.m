%% short version
if (log_letter_status(ii-1)==5 || log_letter_status(ii-1)==10)
    itmp=[5 10 4 7 1];
elseif log_letter_status(ii-1)==4
    itmp=[4 1 5 10 7];
elseif log_letter_status(ii-1)==7
    itmp=[7 1 5 10 4];
elseif (log_letter_status(ii-1)==1)
    itmp=[1 5 10 4 7];
else
    if (log_letter_status(ii-2)==5 || log_letter_status(ii-3)==5 || log_letter_status(ii-2)==10 || log_letter_status(ii-3)==10)
        itmp=[4 7 5 10 1];
    elseif (log_letter_status(ii-2)==7 || log_letter_status(ii-3)==7)
        itmp=[1 5 10 7 4];
    elseif (log_letter_status(ii-2)==4 || log_letter_status(ii-3)==4)
        itmp=[1 5 10 4 7];
    else
        itmp=[5 10 1 4 7];
    end
end
%% Full version
% % % if (log_letter_status(ii-1)==5 || log_letter_status(ii-1)==10)
% % %     itmp=[5 10 4 7 6 1 3];
% % % elseif log_letter_status(ii-1)==4
% % %     itmp=[4 1 5 10 7 6 3];
% % % elseif log_letter_status(ii-1)==7
% % %     itmp=[7 1 5 10 4 6 3];
% % % elseif log_letter_status(ii-1)==6
% % %     itmp=[6 1 5 10 4 7 3];
% % % elseif log_letter_status(ii-1)==3
% % %     itmp=[3 1 5 10 4 7 6];
% % % elseif (log_letter_status(ii-1)==1)
% % %     itmp=[1 5 10 4 7 6 3];
% % % else
% % %     if (log_letter_status(ii-2)==5 || log_letter_status(ii-3)==5 || log_letter_status(ii-2)==10 || log_letter_status(ii-3)==10)
% % %         itmp=[4 7 6 5 10 1 3];
% % %     elseif (log_letter_status(ii-2)==7 || log_letter_status(ii-3)==7)
% % %         itmp=[1 5 10 7 4 6 3];
% % %     elseif (log_letter_status(ii-2)==4 || log_letter_status(ii-3)==4)
% % %         itmp=[1 5 10 4 7 6 3];
% % %     elseif (log_letter_status(ii-2)==6 || log_letter_status(ii-3)==6)
% % %         itmp=[1 5 10 6 4 7 3];
% % %     else
% % %         itmp=[5 10 1 4 7 6 3];
% % %     end
% % % end